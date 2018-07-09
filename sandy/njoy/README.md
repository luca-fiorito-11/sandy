#### <a name="description"></a>Description
The ```sandy.njoy```module effectively generates and runs NJOY input files with only a minimal interaction from the user side, thus making the nuclear data processing accessible to a larger community.

With only few entries provided in the command line, ```sandy.njoy``` automatically determines the best NJOY input options for a succesfull run.
On top of that, several file-related parameters --- e.g. MAT number---  are read and imported from the desired evaluated file without having the user dig into the details of a very human-unfriendly format such as ENDF-6.

- [Before starting](#before_staring)
- [NJOY versions](#versions)
- [Usage](#usage)
- [Producing PENDF files](#PENDF)
- [Producing GENDF files](#PENDF)
- [Producing ERRORR files](#ERRORR)
- [Producing ACE files](#ACE)


#### <a name="before_starting"></a> Before starting

```sandy.njoy``` looks  for the executable file `njoy2016` in the ```PATH``` environmetal variable. If you want to use a different NJOY executable, then you must specify the desired file using the dedicated command line argument.

#### <a name="versions"></a> NJOY versions

Notice that ```sandy.njoy``` only supports [NJOY-2016](https://njoy.github.io/NJOY2016/).

#### <a name="usage"></a> Usage

```sandy.njoy``` is run as a ```python```script using the ```-m``` command line option

```bash
$ python -m sandy.njoy
```

For an overview of the ```sandy.njoy``` usage and acceptable arguments, type

```bash
$ python -m sandy.njoy --help
```
```
usage: python -m sandy.njoy [-h] [-P PENDF] [-G GENDF] [-p] [-a] [-g] [-e]
                            [-H] [--temps TEMPS [TEMPS ...]]
                            [--sig0 SIG0 [SIG0 ...]]
                            [--kerma [KERMA [KERMA ...]]] [--err ERR]
                            [--free-gas] [--ptable] [--gaspr] [--ign IGN]
                            [--iwt IWT] [--igne IGNE] [--iwte IWTE]
                            [--suffixes .XX [.XX ...]] [-V VERBOSE] [-v]
                            tape

positional arguments:
  tape                  ENDF-6 format file

optional arguments:
  -h, --help            show this help message and exit
  -P PENDF, --pendftape PENDF
                        processed PENDF format file
  -G GENDF, --gendftape GENDF
                        processed GENDF format file
  -p, --pendf           produce PENDF file
  -a, --ace             produce ACE file
  -g, --gendf           produce GENDF file
  -e, --errorr          produce ERRORR file
  -H, --hendf           produce HENDF file
  --temps TEMPS [TEMPS ...]
                        list of temperature values
  --sig0 SIG0 [SIG0 ...]
                        list of dilution values
  --kerma [KERMA [KERMA ...]]
                        list of partial kermas (default=[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447])
  --err ERR             fractional tolerance for RECONR and BROADR (default=0.001)
  --free-gas            compute thermal cross section for free-gas scatterer (THERMR)
  --ptable              compute probability tables (PURR)
  --gaspr               compute gas-production cross sections (GASPR)
  --ign IGN             neutron group structure option for GROUPR
                         - file : read from file
                         - 2 : csewg 239-group structure [default]
                         - 3 : lanl 30-group structure
                         - 4 : anl 27-group structure
                         - 5 : rrd 50-group structure
                         - 6 : gam-i 68-group structure
                         - 7 : gam-ii 100-group structure
                         - 8 : laser-thermos 35-group structure
                         - 9 : epri-cpm 69-group structure
                         - 10 : lanl 187-group structure
                         - 11 : lanl 70-group structure
                         - 12 : sand-ii 620-group structure
                         - 13 : lanl 80-group structure
                         - 14 : eurlib 100-group structure
                         - 15 : sand-iia 640-group structure
                         - 16 : vitamin-e 174-group structure
                         - 17 : vitamin-j 175-group structure
                         - 19 : ecco 33-group structure
                         - 20 : ecco 1968-group structure
                         - 21 : tripoli 315-group structure
                         - 22 : xmas lwpc 172-group structure
                         - 23 : vit-j lwpc 175-group structure
                        predefined group structures:
                         - scale_238
  --iwt IWT             Weight function option for GROUPR
                         - file : read from file
                         - 2 : constant
                         - 3 : 1/e
                         - 5 : epri-cell lwr
                         - 6 : (thermal) -- (1/e) -- (fission + fusion)
                         - 7 : same with t-dep thermal part
                         - 8 : thermal--1/e--fast reactor--fission + fusion
                         - 9 : claw weight function
                         - 10 : claw with t-dependent thermal part
                         - 11 : vitamin-e weight function (ornl-5505)
                         - 12 : vit-e with t-dep thermal part
                        predefined functions:
                         - jaea_fns_175
  --igne IGNE           neutron group structure option for ERRORR (same as --ign)
  --iwte IWTE           weight function option for ERRORR (same as --iwt)
  --suffixes .XX [.XX ...]
                        suffixes for ACE files, as many as temperature values (default = None).
  -V VERBOSE, --verbose VERBOSE
                        set verbosity level (default = 1)
  -v, --version
```

#### <a name="PENDF"></a> Producing PENDF files
Let's take a look at the various ```sandy.njoy``` options to produce a pointwise ENDF (PENDF) file.

##### CROSS SECTION RECONSTRUCTION
First off, a PENDF file is produced by running ```sandy.njoy``` with option ```-p``` or ```--pendf```.
Cross sections are reconstructed and linearized using the RECONR module of NJOY.
```bash
$ python -m sandy.njoy ${file} -p
```
The PENDF output will be named ```${file}.pendf```.

##### DOPPLER-BROADENING
Cross sections are Doppler-broadened at one or more temperatures using the BROADR module of NJOY.
For example, you can include Doppler-broadening for, say, two temperatures (300 and 600 K), using option ```--temps``` like this:
```bash
$ python -m sandy.njoy ${file} -p \
--temps 300 600
```

##### URR TREATMENT
To produce effective self-shielded cross sections for resonance reactions in the unresolved energy range --- using the UNRESR module of NJOY --- you must define the desired levels of dilutions with keyword ```--sig0```, as
```bash
$ python -m sandy.njoy ${file} -p \
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0
```

##### HEATING / DAMAGE ENERGY
Nuclear heating and/or damage energy can be computed with the HEATR module of NJOY by using keyword ```--kerma```.
```bash
$ python -m sandy.njoy ${file} -p \
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma
```
The ```--kerma``` instruction accepts a list of heating and damage MT reaction numbers that will be computed and added to the PENDF file.
If ```--kerma``` is provided without additional entries, the default ```--kerma 302 303 304 318 402 442 443 444 445 446 447``` is used.

The allowed MT reaction numbers with their meaning is reported in the table below.
| MT | Description |
| ------ | ------ |
| 301 | KERMA total (energy balance) |
| 302 | KERMA elastic |
| 303 | KERMA non-elastic |
| 304 | KERMA inelastic |
| 318 | KERMA fission  |
| 401 | KERMA disappearance  |
| 402 | KERMA radiative capture |
| 403 | KERMA proton production  |
| 407 | KERMA alpha production |
| 442 | total photon KERMA contribution |
| 443 | total kinematic KERMA (kinematic Limit)  |
| 444 | total damage energy production cross section  |
| 445 | elastic damage energy production cross section  |
| 446 | inelastic damage energy production cross section |
| 447 | neutron disappearance damage energy production cross section |

##### THERMAL SCATTERING
Thermal cross sections for a free-gas scatterer can be computed by the THERMR module of NJOY by adding the ```--free-gas``` instruction.
```bash
$ python -m sandy.njoy ${file} -p \
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma \
--free-gas
```
The inelastic reaction is added to the PENDF file with MT reaction number 221.

##### PROBABILITY TABLES
Option ```--ptable``` adds probability tables --- calculated with the module PURR of NJOY --- to the PENDF file.
```bash
$ python -m sandy.njoy ${file} -p \
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma \
--free-gas \
--ptable
```

##### GAS PRODUCTION
Gas production cross sections can be added to the PENDF file --- with the module GASPR of NJOY --- using instruction ```--gaspr```.
```bash
$ python -m sandy.njoy ${file} -p \
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma \
--free-gas \
--ptable \
--gaspr
```

#### <a name="GENDF"></a> Producing GENDF files
GENDF or *groupwise-ENDF* files are processed files produced by the GROUPR module of NJOY. Such files are written according to a modified ENDF-6 format conceived for storing group average nuclear data.
To produce GENDF files by running the ```sandy.njoy``` module, you can try and use instruction ```-g``` or ```--gendf```, like this:
```bash
$ python -m sandy.njoy ${file} -g
```

In this case, an error messages wil be raised stating that a PENDF file is missing (NJOY needs a ENDF-6 and a PENDF file to produce a GENDF output).
You can provide the PENDF file either usign keyword ```--pendftape ${pendf_file}```, as
```bash
$ python -m sandy.njoy ${file} -g \
--pendftape ${pendf_file}
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0
```

or including ```-g``` to the PENDF production sequence:
```bash
$ python -m sandy.njoy ${file} -pg \ # process both pendf and gendf
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma \
--free-gas \
--ptable \
--gaspr
```
The generated GENDF file will be called ```${file}.gendf```.

Notice that, independently on the way the GENDF file is produced, ```sandy.njoy``` requires an explicit specification of the temperature and dilution values.

##### GROUP STRUCTURE
As an extra instruction, you can define a neutron group structure for the output GENDF file other than the default CSWEG 239-group structure.
Then, add argument ```--ign``` with the desired structure number according to the NJOY manual, e.g.
```bash
$ python -m sandy.njoy ${file} -pg \ # process both pendf and gendf
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma \
--free-gas \
--ptable \
--gaspr \
--ign 19      # produce GENDF output file on the ECCO 33-group structure
```

A tabulated group structure can be read from a text file if the filename replaces the NJOY structure number, ```--ign ${text_filename}```.
In such a case, each group boundary should be provided in eV as a new line of the file, like this:
```ascii
1E-5
1E-3
0.1
1.0000
100
1e4
1e7
```

Alternatively, precompiled structures that are not implemented in NJOY can be used,
For example, to use the 238-group structure included in the SCALE suite, just type ```--ign scale_238```.

##### WEIGHTING FUNCTION
Analogously, the default weighting function (thermal -- 1/E -- fission + fusion) specified for the group collapsing can be modified using keyword ```--iwt``` followed by one of the dedicated numbers defined in the NJOY manual, for instance
```bash
$ python -m sandy.njoy ${file} -pg \ # process both pendf and gendf
--temps 300 600 \
--sig0 1E10 1E4 1E3 1E2 1E1 1E0 \
--kerma \
--free-gas \
--ptable \
--gaspr \
--ign 19 \
--iwt 2     # use a constant function
```

The weighting function can be also read from a text file typing ```--iwt ${text_filename}```.
Then, the requested format for the file is the following:
 - two columns, the first for the energy (in eV) and the second for the weight;
 - the weight in line *i* is considered constant between the two energies in lines *i* and *i+1*.

```ascii
1E-5  1
1E-3  3
0.1   3.5
1.0000  4
100  4.5
1e4  5
1e7  0
```

Precompiled data also exist for weighting functions.

#### <a name="ERRORR"></a> Producing ERRORR files
The ERRORR module of NJOY produces multigroup covariance matrices by convoluting the covariance information available in ENDF-6 evaluated files and the corresponding central values.
To produce ERRORR files with ```sandy.njoy```, add option ```-e``` or ```--errorr``` to the PENDF or GENDF sequence:
```bash
$ python -m sandy.njoy ${file} -pe \ # process pendf and errorr
--temps 300 600
```
The generated ERRORR file will be called ```${file}.errorr```.

Alternatively, run the module providing the PENDF file as an extra argument:
```bash
$ python -m sandy.njoy ${file} -e \ # process errorr
--temps 300 600 \
--pendftape ${pendf_file}
```

For both cases, the temperature values must be given explicitely using keyword ```--temps```.

##### GROUP STRUCTURE
As for the production of GENDF files, a user-defined group structure can be provided using keyword ```--igne```.
For such a keyword the same rules as or ```--ign``` apply.

##### WEIGHTING FUNCTION
You can also supply a weighting funcion using command ```--iwte```, for which the rules defined for ```--iwt``` apply.

#### <a name="ACE"></a> Producing ACE files

ACE (A Compact ENDF) formatted files are used by many continous-energy Monte Carlo transport codes including MCNP.
You can process ACE files running ```sandy.njoy```  with option ```-a``` or ```--ace```, as
```bash
$ python -m sandy.njoy ${file} -a \ # process ace
--temps 300 600 \
--pendftape ${pendf_file}
```
The generated ACE file will be called ```${file}.ace```.

If you do not have a suitable PENDF file available, simply reprocess it with ```sandy.njoy``` and attach option ```-a``` to the sequence.
```bash
$ python -m sandy.njoy ${file} -pa \ # process pendf and ace
--temps 300 600
```

Remember that the temperature specification is mandatory, as well as the respective suffixes used in the xsdir file.
Suffixes can be specified as it follows:
```bash
$ python -m sandy.njoy ${file} -pa \ # process pendf and ace
--temps 300 600 \
--suffixes .03 .06
```
Keyword ```--suffixes``` takes as many entries as given with ```--temps```.