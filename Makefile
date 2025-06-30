#Per runnare il docker; fai exec dopo che il docker è start

DIR=C:\Users\tiziana.liberati\dockers\ProgettoTesi

image:
	docker build . -t kemopenfhe

run:
	docker container rm kemopenfhe & docker run --name kemopenfhe -v "$(DIR)/src:/app/target" -v "$(DIR)/openfhe-development:/usr/local/include/openfhe" -it kemopenfhe bash

exec:
	docker start kemopenfhe & docker exec -it kemopenfhe bash