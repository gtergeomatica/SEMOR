### Setup SEMOR

### Compile
Commands to be executed inside the SEMOR folder:
```
  mkdir build
  cd build
  cmake ../
  make
  cd ..
  cd RTKLIB-b34e/app/consapp/rtkrcv/gcc/
  make
  cd ../../../../..
  cd RTKLIB-b34e/app/consapp/str2str/gcc/
  make
```
### Execute
Commands to be executed inside the SEMOR folder:
```
  cd bin
  ./semor
```

### Configuration
After the first execution of SEMOR a default semor.conf is created.
Update it based on your needs.
  
In order to stop SEMOR you need to send 'q' to the terminal.
  
The output is written in the output.txt file inside the root folder (SEMOR).

The output is not written immediately, SEMOR have to wait for input from rtkrcv instances.
