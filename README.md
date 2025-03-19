# LTTC Homework CIS and TDHF

This project contains the electronic structure theory program for calculating the restricted Hartree-Fock energy as well as the excitation energies by means of CIS and TDHF methods.

## Project Structure

This project is composed of several folders. For instance, the `examples` folder has the example files, the `input` folder contains the information about the basis sets and other relevant data for simulations, the `include` folder has a file with the important parameters, the `report` folder contains the pdf and TeX files of the report, the `src` folder has the source code, the `test` folder has the example output as well as the `LICENSE` document.

```
└── 📁LTTC-Homework--CIS-TDHF
    └── 📁examples
    └── 📁include
        └── parameters.h
    └── 📁input
        └── basis
        └── molecule
        └── options
    └── 📁int
        └── ERI.dat
        └── Kin.dat
        └── Nuc.dat
        └── Ov.dat
    └── 📁report
        └── report.pdf
        └── report.tex
    └── 📁src
        └── antisymmetrize_ERI.f90
        └── AO_to_MO.f90
        └── CIS.f90
        └── EST.f90
        └── Makefile
        └── 📁obj
        └── orthogonalization_matrix.f90
        └── print_RHF.f90
        └── read_basis.f90
        └── read_geometry.f90
        └── read_integrals.f90
        └── read_molecule.f90
        └── RHF.f90
        └── spatial_to_spin_ERI.f90
        └── spatial_to_spin_MO_energy.f90
        └── TDHF.f90
        └── utils.f90
        └── wrap_lapack.f90
    └── 📁test
        └── Be_cc-pvdz.out
        └── Be_cc-pvtz.out
        └── He_cc-pvdz.out
        └── He_cc-pvtz.out
        └── Ne_cc-pvdz.out
        └── Ne_cc-pvtz.out
    └── GoEST
    └── Instructions.pdf
    └── LICENSE
    └── README.md
```

## Required software

Ensure you have the following installed on your system:
- `gfortran` (version 14.2.0 was tested)
- `make` (version 3.81 was tested)

## Installation

1. First, clone the repository:
    ```sh
    git clone https://github.com/almakhmudov/LTTC-Homework--CIS-TDHF.git
    cd LTTC-Homework--CIS-TDHF
    ```

2. Compile the program using `make`:
    ```sh
    cd src
    make
    ```

3. In order to delete the executable run the following:
    ```sh
    cd src
    make clean
    ```

## Test

The program was tested on MacOS Sequoia 15.2. To test whether the installation was successful, you could run the calculations for He, Ne, and Be using the cc-pvdz and cc-pvtz basis sets and compare the obtained results to the expected output in the `test` folder. Alternatively, one could check the obtained values with the table below.

<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  overflow:hidden;padding:10px 5px;word-break:normal;}
.tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
.tg .tg-c3ow{border-color:inherit;text-align:center;vertical-align:top}
.tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}
</style>
<table class="tg"><thead>
  <tr>
    <th class="tg-0pky"></th>
    <th class="tg-c3ow" colspan="2">CIS</th>
    <th class="tg-c3ow" colspan="2">TDHF</th>
  </tr></thead>
<tbody>
  <tr>
    <td class="tg-0pky"></td>
    <td class="tg-c3ow">cc-pvdz</td>
    <td class="tg-c3ow">cc-pvtz</td>
    <td class="tg-c3ow">cc-pvdz</td>
    <td class="tg-c3ow">cc-pvtz</td>
  </tr>
  <tr>
    <td class="tg-c3ow">&nbsp;&nbsp;He  </td>
    <td class="tg-c3ow">51.947</td>
    <td class="tg-c3ow">31.807</td>
    <td class="tg-c3ow">51.577</td>
    <td class="tg-c3ow">31.608</td>
  </tr>
  <tr>
    <td class="tg-c3ow">&nbsp;&nbsp;Ne  </td>
    <td class="tg-c3ow">49.010</td>
    <td class="tg-c3ow">33.210</td>
    <td class="tg-c3ow">48.871</td>
    <td class="tg-c3ow">33.174</td>
  </tr>
  <tr>
    <td class="tg-c3ow">&nbsp;&nbsp;Be  </td>
    <td class="tg-c3ow">5.295</td>
    <td class="tg-c3ow">5.143</td>
    <td class="tg-c3ow">4.993</td>
    <td class="tg-c3ow">4.867</td>
  </tr>
</tbody>
</table>

## Usage

In order to run the program, one should first change into the main directory, namely `LTTC-Homework--CIS-TDHF`. From there one could execute the example line as mentioned below. Following atoms and molecules are supported: He, Ne, Be, H2O together with the cc-pvdz and cc-pvtz basis sets.

Example:
```sh
./GoEST <atom or molecule> <basis set>
./GoEST He cc-pvdz
```

> [!IMPORTANT]
> This program works only with the aforementioned atoms, molecules and basis sets combinations. Otherwise it won't produce any results.

## Acknowledgments
This project is based upon data and instructions provided by Pina Romaniello.