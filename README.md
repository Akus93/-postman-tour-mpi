**[WE][Inf][mgr][PrP]**
# **Problem chińskiego listonosza**


## Nowy sposób kompilacji i uruchamiania.
!!! Nie testowany w PCSSIE !!! W razie problemów można ręcznie (po staremu) kompilować/uruchamiać.

#### Kompilacja
>bash ./compile.sh

#### Uruchomienie
> bash ./start.sh {opcjonalnie liczba procesów, domyślnie 10}

**W skryptach nie ma instrukcji odpowiedzialnych za przygotowanie srodowiska w PCSSie więc trzeba albo ręcznie albo dopisać do skryptu**

srun --pty /bin/bash

module load openmpi

nano mpi_hello_world.c

mpiCC mpi_hello_world.c -o mpi_hello_world.mpi

salloc -N 10 mpirun mpi_hello_world.mpi



Z powodu korzystania z nowszych wersji języka możliwe że przy ręcznej kompilacji wymagane będzie dodanie parametru --std=c++11
