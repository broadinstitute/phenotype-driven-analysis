import tempfile

import gspread
import os
import sqlite3
import subprocess
import numpy as np
import pysam
from google.oauth2.service_account import Credentials
from google.cloud import storage
from pprint import pformat

import hail as hl
#import hailtop.fs as hfs

"""
value_render_option: can be "FORMATTED_VALUE", "UNFORMATTED_VALUE", or "FORMULA":

    FORMATTED_VALUE: For example, if A1 is 1.23 and A2 is =A1 and formatted as currency, then A2 would return "$1.23".
    UNFORMATTED_VALUE: For example, if A1 is 1.23 and A2 is =A1 and formatted as currency, then A2 would return the number 1.23.
    FORMULA: The reply will include the formulas. For example, if A1 is 1.23 and A2 is =A1 and formatted as currency, then A2 would return "=A1".
     
https://developers.google.com/sheets/api/reference/rest/v4/ValueRenderOption
"""

VALUE_RENDER_OPTION__FORMATTED_VALUE = "FORMATTED_VALUE"
VALUE_RENDER_OPTION__UNFORMATTED_VALUE = "UNFORMATTED_VALUE"
VALUE_RENDER_OPTION__FORMULA = "FORMULA"

# spreadsheets must be shared with 733952080251-compute@developer.gserviceaccount.com
_GOOGLE_CREDETIALS_JSON_PATH = os.path.expanduser('~/.config/gcloud/seqr-project-0cb2b89f436f.json')
_GSPREAD_CLIENT = None

def get_storage_client():
    #  based on https://googleapis.dev/python/google-api-core/latest/auth.html
    return storage.Client() # use credentials from gcloud auth application-default login


def get_spreadsheet(spreadsheet_name):
    global _GSPREAD_CLIENT
    if _GSPREAD_CLIENT is None:
        creds = Credentials.from_service_account_file(
            _GOOGLE_CREDETIALS_JSON_PATH,
            scopes=[
                'https://www.googleapis.com/auth/spreadsheets',
                'https://www.googleapis.com/auth/drive.file',
                'https://www.googleapis.com/auth/drive',
            ]
        )

        _GSPREAD_CLIENT = gspread.authorize(creds)

    spreadsheet = _GSPREAD_CLIENT.open(spreadsheet_name)

    return spreadsheet


_bam_header_date_cache_db = None
def _connect_to_bam_header_date_cache_db():
    global _bam_header_date_cache_db
    if _bam_header_date_cache_db is not None:
        return _bam_header_date_cache_db

    cache_db_path = '~/code/sample_metadata/metadata/bam_header_date_cache.db'
    _bam_header_date_cache_db = sqlite3.connect(
        os.path.expanduser(cache_db_path),
        isolation_level=None,
        cached_statements=0)
    print("Connected to cache_db: ", cache_db_path)
    try:
        _bam_header_date_cache_db.execute("CREATE TABLE cache (path, bam_header_date)").close()
        _bam_header_date_cache_db.execute("CREATE UNIQUE INDEX cache_index ON cache (path)").close()

    except sqlite3.OperationalError as e:
        if "already exists" not in str(e):
            print("ERROR:", e)

    return _bam_header_date_cache_db


def get_date_from_bam_header(bam_or_cram_path, gcloud_project="bw2-rare-disease"):
    # check the cache
    cache_db = _connect_to_bam_header_date_cache_db()

    cursor = cache_db.execute(f"SELECT path, bam_header_date FROM cache where path=?", (bam_or_cram_path,))
    try:
        cached_bam_path, bam_header_date = next(cursor)
        return bam_header_date
    except StopIteration:
        pass
    finally:
        cursor.close()

    output = subprocess.check_output(f"gsutil -u {gcloud_project} cat %s | samtools view -H - | grep ^@RG | head -n 1" % bam_or_cram_path, shell=True, encoding="UTF-8", stderr=subprocess.DEVNULL)
    read_group_annotations = {}
    for i, rg_field in enumerate(output.rstrip().split("\t")):
        if i == 0:
            continue  # skip @RG prefix
        key = rg_field[:rg_field.find(":")]
        value = rg_field[rg_field.find(":")+1:]
        read_group_annotations[key] = value

    try:
        bam_header_date = read_group_annotations['DT'][:7]
    except Exception as e:
        message = f"Unable to find date ('DT') key in read_group_annotations: {pformat(read_group_annotations)}: {e}"
        print(f"ERROR: {message}")
        return ""
        #raise ValueError(message)

    # update cache
    #print(f"Adding {bam_or_cram_path} => {genome_version} to sqlite cache")
    if bam_header_date:
        cache_db.execute(f"INSERT INTO cache VALUES (?, ?)", (bam_or_cram_path, bam_header_date)).close()

    return bam_header_date


_cram_and_bam_path_genome_version_cache_db = None
def _connect_to_cram_and_bam_path_genome_version_cache_db():
    global _cram_and_bam_path_genome_version_cache_db
    if _cram_and_bam_path_genome_version_cache_db is not None:
        return _cram_and_bam_path_genome_version_cache_db

    cache_db_path = 'cram_and_bam_path_genome_version_cache.db'
    _cram_and_bam_path_genome_version_cache_db = sqlite3.connect(
        os.path.expanduser(cache_db_path),
        isolation_level=None,
        cached_statements=0)
    print("Connected to cache_db: ", cache_db_path)
    try:
        _cram_and_bam_path_genome_version_cache_db.execute("CREATE TABLE cache (path, genome_version INT(8))").close()
        _cram_and_bam_path_genome_version_cache_db.execute("CREATE UNIQUE INDEX cache_index ON cache (path)").close()

    except sqlite3.OperationalError as e:
        if "already exists" not in str(e):
            print("ERROR:", e)

    return _cram_and_bam_path_genome_version_cache_db


def get_genome_version_from_bam_or_cram_header(bam_or_cram_path, gcloud_project="bw2-rare-disease"):
    # check the cache
    cache_db = _connect_to_cram_and_bam_path_genome_version_cache_db()

    cursor = cache_db.execute(f"SELECT path, genome_version FROM cache where path=?", (bam_or_cram_path,))
    try:
        cached_bam_path, cached_genome_version = next(cursor)
        return cached_genome_version
    except StopIteration:
        pass
    finally:
        cursor.close()

    # get genome version from file header
    output = subprocess.check_output(
        f"gsutil -u {gcloud_project} cat %s | samtools view -H - | grep @SQ | head -n 3" % bam_or_cram_path, shell=True, encoding="UTF-8", stderr=subprocess.DEVNULL)
    genome_version = None
    if "GRCh37" in output or "Homo_sapiens_assembly19.fasta" in output:
        genome_version = 37
    elif "GRCh38" in output or "Homo_sapiens_assembly38.fasta" in output:
        genome_version = 38
    else:
        print(f"WARNING: Unable to parse genome version from header lines in {bam_or_cram_path}: {output}")

    # update cache
    #print(f"Adding {bam_or_cram_path} => {genome_version} to sqlite cache")
    if genome_version:
        print(f"Got genome version hg{genome_version} from {bam_or_cram_path} header")
        cache_db.execute(f"INSERT INTO cache VALUES (?, ?)", (bam_or_cram_path, genome_version)).close()

    return genome_version


_cram_and_bam_path_fragment_length_stats_cache_db = None
def _connect_to_cram_and_bam_path_fragment_length_stats_cache_db():
    global _cram_and_bam_path_fragment_length_stats_cache_db
    if _cram_and_bam_path_fragment_length_stats_cache_db is not None:
        return _cram_and_bam_path_fragment_length_stats_cache_db

    cache_db_path = 'cram_and_bam_path_fragment_length_stats_cache.db'
    _cram_and_bam_path_fragment_length_stats_cache_db = sqlite3.connect(
        os.path.expanduser(cache_db_path),
        isolation_level=None,
        cached_statements=0)
    print("Connected to cache_db: ", cache_db_path)
    try:
        _cram_and_bam_path_fragment_length_stats_cache_db.execute("CREATE TABLE cache (path, read_length FLOAT, coverage FLOAT, fragment_length_mean FLOAT, fragment_length_stdev FLOAT)").close()
        _cram_and_bam_path_fragment_length_stats_cache_db.execute("CREATE UNIQUE INDEX cache_index ON cache (path)").close()

    except sqlite3.OperationalError as e:
        if "already exists" not in str(e):
            print("ERROR:", e)

    return _cram_and_bam_path_fragment_length_stats_cache_db


def get_fragment_length_stats_from_bam_or_cram(bam_or_cram_path, reference_fasta_path, gcloud_project="bw2-rare-disease", n_reads_to_process=1_000_000, force=False):
    # check the cache
    cache_db = _connect_to_cram_and_bam_path_fragment_length_stats_cache_db()

    cached_bam_path = None
    cursor = cache_db.execute(f"SELECT path, read_length, coverage, fragment_length_mean, fragment_length_stdev FROM cache where path=?", (bam_or_cram_path,))
    try:
        cached_bam_path, read_length, coverage, fragment_length_mean, fragment_length_stdev = next(cursor)
        if not force:
            return read_length, coverage, fragment_length_mean, fragment_length_stdev
    except StopIteration:
        pass
    finally:
        cursor.close()

    temp_bam = tempfile.NamedTemporaryFile("wb", suffix=".bam", delete=True)

    # compute read_length, coverage, fragment_length_mean, fragment_length_stdev from the bam file
    print(f"Downloading {n_reads_to_process} reads from {bam_or_cram_path} to {temp_bam.name}")
    subprocess.check_call(
            f"gsutil -u {gcloud_project} cat {bam_or_cram_path} | "
            f"samtools view -h - | "
            f"head -n {n_reads_to_process+2000} | "  # add 2000 to take header into account
            f"samtools view -b - > {temp_bam.name} && "
            f"samtools index {temp_bam.name}",
        shell=True, stderr=subprocess.DEVNULL)

    read_length = 0
    total_bases = 0
    start_coord = 10**9
    end_coord = 0
    fragment_lengths = []
    for r in pysam.AlignmentFile(temp_bam.name, reference_filename=reference_fasta_path):
        if r.is_unmapped or r.is_secondary or not r.is_proper_pair or r.mapping_quality < 50:
            continue
        if r.reference_name != r.next_reference_name:
            continue
        fragment_length = abs(r.next_reference_start - r.pos + 1)
        if 10 < fragment_length < 1500:
            start_coord = min(start_coord, r.pos)
            end_coord = max(end_coord, r.pos + read_length)
            read_length = max(read_length, r.query_length)
            total_bases += r.query_length
            fragment_lengths.append(fragment_length)
    temp_bam.close()

    coverage = total_bases/(end_coord-start_coord) if end_coord-start_coord > 0 else 0
    fragment_length_mean = np.mean(fragment_lengths)
    fragment_length_stdev = np.std(fragment_lengths)

    # update cache
    print(f"Adding {bam_or_cram_path} => read_length: {read_length}, coverage: {coverage}, mean: {fragment_length_mean}, stdev: {fragment_length_stdev} to sqlite cache")
    if cached_bam_path:
        cache_db.execute(f"UPDATE cache SET read_length=?, coverage=?, fragment_length_mean=?, fragment_length_stdev=? WHERE path=?", (read_length, coverage, fragment_length_mean, fragment_length_stdev, bam_or_cram_path)).close()
    else:
        cache_db.execute(f"INSERT INTO cache VALUES (?, ?, ?, ?, ?)", (bam_or_cram_path, read_length, coverage, fragment_length_mean, fragment_length_stdev)).close()

    return read_length, coverage, fragment_length_mean, fragment_length_stdev


_hadoop_cache_db = None

def _connect_to_hadoop_cache_db():
    global _hadoop_cache_db
    if _hadoop_cache_db is not None:
        return _hadoop_cache_db

    cache_db_path = 'hadoop_cache.db'
    _hadoop_cache_db = sqlite3.connect(
        os.path.expanduser(cache_db_path),
        isolation_level=None,
        cached_statements=0)
    print("Connected to cache_db: ", cache_db_path)
    try:
        _hadoop_cache_db.execute("CREATE TABLE cache (path, operation, result)").close()
        _hadoop_cache_db.execute("CREATE UNIQUE INDEX cache_index ON cache (path)").close()
    except sqlite3.OperationalError as e:
        if "already exists" not in str(e):
            print("ERROR:", e)

    return _hadoop_cache_db


def hadoop_exists_using_cache(path, double_check_if_cache_says_yes=False, double_check_if_cache_says_no=False):
    # check the cache
    cache_db = _connect_to_hadoop_cache_db()
    operation = "exists"
    cursor = cache_db.execute(f"SELECT path, result FROM cache where path=? and operation=?", (path, operation,))
    file_exists = None
    try:
        path, result = next(cursor)
        if result not in ("True", "False"):
            raise ValueError(f"Unexpected value in  cache for {path} {operation}: {result}")
        file_exists = result == "True"
        if file_exists and not double_check_if_cache_says_yes:
            return True
        if not file_exists and not double_check_if_cache_says_no:
            return False

    except StopIteration:
        pass
    finally:
        cursor.close()

    # do operation
    result = hl.hadoop_exists(path)
    #result = hfs.exists(path)

    # update cache
    #print(f"Adding {bam_or_cram_path} => {genome_version} to sqlite cache")
    if file_exists is not None:
        if str(file_exists) != str(result):
            print(f"Updating {path} from {file_exists} to {result}")
        cache_db.execute(f"UPDATE cache SET result=? WHERE path=? AND operation=?", (str(result), path, operation)).close()
    else:
        cache_db.execute(f"INSERT INTO cache VALUES (?, ?, ?)", (path, operation, str(result))).close()

    return result


def hadoop_file_size_using_cache(path, force=False):
    """Returns the file size in bytes"""

    cached_file_size = None

    # check the cache
    cache_db = _connect_to_hadoop_cache_db()
    operation = "size"
    cursor = cache_db.execute(f"SELECT path, result FROM cache where path=? and operation=?", (path, operation,))
    try:

        path, cached_file_size = next(cursor)
        try:
            cached_file_size = int(cached_file_size)
        except:
            raise ValueError(f"Unable to convert file size to integer {path} {operation}: {cached_file_size}")

        if not force:
            return cached_file_size

    except StopIteration:
        pass
    finally:
        cursor.close()

    # do operation
    file_stat = hl.hadoop_stat(path)
    file_size = int(file_stat['size_bytes'])

    # update cache
    print(f"Got file size: {file_size:,d} bytes for {path}")
    if cached_file_size is not None:
        cache_db.execute(f"UPDATE cache SET result=? WHERE path=? AND operation=?", (str(file_size), path, operation)).close()
    else:
        cache_db.execute(f"INSERT INTO cache VALUES (?, ?, ?)", (path, operation, str(file_size))).close()

    return file_size

