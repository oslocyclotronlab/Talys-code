# Talys-code
This repository is created to store the source code of the TALYS package. Since the creators of the package does not provide any bug-tracker, change history, etc., we will have a private repository to keep track of changes here. If you recive a 'bug-fix' or similar, please create a new branch and add the change there. That way we might be able to keep track of what is fixed, what has been changed etc. Also, we should use the bug-tracker when anyone finds a bug. Currently the default branch is the 'newest' 1.9 version with one or more bugfixes. Please add earlier versions such as 1.8 in new branches.

## Notes:
I've added some changes in this branch:
* Improved build system - I've added a CMake script to build the package. This will make it possible to make a simple Homebrew formula
* Added the ability to give density tables as input without having to edit the database

### Density table
To use tabulated level densities from a file without modifying the database you can use the following input:
```
densfile Z A filename
```
where Z is the atomic number, A the mass number of the nucleus' density to replace. `filename` is the path to the density table. The table should be formatted as in the tabulated level density files in the database. If you want to give separate tables for each parity add the following line to the input file:
```
ldmodel Z A k
```
where `k` is either 5 or 6.

An example input file with density file are shown in the `test` folder.

#### Note: There may be some issues with the parity options. Needs more investigation.

### Note: I have not included the database part of the package as that is a huge package of data.
