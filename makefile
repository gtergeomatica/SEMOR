default: semor

semor.o:
	gcc -c semor.c -g -o bin/semor.o

client.o:
	gcc -c src/client.c -g -o bin/client.o

imu_log.o:
	gcc -c src/imu_log.c -g -o bin/imu_log.o

semor: semor.o client.o imu_log.o
	gcc bin/semor.o bin/client.o bin/imu_log.o -o bin/semor -lm -lpthread

clean:
	-rm -f bin/semor.o
	-rm -f bin/semor
	-rm -f bin/client.o