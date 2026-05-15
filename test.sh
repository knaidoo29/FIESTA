coverage erase
coverage run --parallel-mode --source=fiesta -m pytest tests/normal_tests
mpirun -n 3 coverage run --parallel-mode --source=fiesta -m pytest tests/mpi_tests
coverage combine
coverage report -m