TAG1 = weisburd/exomiser:v1
TAG2 = weisburd/exomiser:latest

all: build push update_sha

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)

update_sha:
	docker pull $(TAG2) 2>&1 | grep Digest | cut -c 9- > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/@sha256:[^\"]*/"@`cat sha256.txt`"/"  ../*.py
