# pyMCD
python code for Multi-Coordinate Driving (MCD) method

## Table of Contents

- [Environmental](#environmental)
- [Input](#input)
- [Output](#output)

## Environmental

- python>=3.7
- [cclib](https://github.com/cclib/cclib)=1.7.1
- [numpy](https://github.com/numpy/numpy)
- [RDKit](https://www.rdkit.org/docs/Install.html)>=2020.09.3
- [scipy](https://github.com/scipy/scipy)

## Input

### Experimental directory Structure
```
example/diels
│
├── R.com       # Geometry of interested molecule
│
├── bond_info   # target bond information
│
└── qc_input    # Input file for Quantum Chemistry package
```
#### 1. R.com
Geometry of target molecule with `charge` and `spin multiplicity`
```
0 1
C 0.00000000 0.00000000 0.00000000
C 0.03101500 0.00321900 1.33695700
C 1.25868100 0.00290300 2.13867800
C 2.49534700 -0.00205700 1.62972600
C 2.26606500 -2.42640100 0.45058000
C 1.15192300 -2.42622200 -0.27662600
H 0.17926100 -2.67118800 0.12875300
H 1.12996800 -2.20777400 -1.33507000
H 2.28637300 -2.67106300 1.50421500
H 3.24383300 -2.20776500 0.04480200
H 3.38321900 -0.04383300 2.24341300
H 2.71059800 0.07480700 0.57265500
H 1.10020000 -0.03373600 3.21984800
H -0.89519400 -0.03471900 1.91664600
H -0.91872700 -0.04169400 -0.56646200
H 0.88142600 0.07857400 -0.62177800
```
#### 2. bond_info
해줘요 ㅋㅋ
```
1 6 1.54 6
4 5 1.54 6
```
#### 3. qc_input
Information containing calculation theory/level, basis set, solvation and other detail options.

**For gaussian,**
```
#N pm6 scf(xqc)
```
**For Orca,**
```
! b3lyp 6-31g(d)
```

### Execute pyMCD
```
python run.py \
	-id <input directory> \
	-sd <save directory> \
	-wd <working directory> \
	-num_relaxation <num_relaxation> \
	-scale <scale> \
	-c <calculator> \
	-u <unit>
```
- `input directory` : The path of directory that input files(e.g. R.com, bond_info, and qc_input) is stored
- `save directory` : The path of directory that output files will be saved
- `working directory` : The path of directory that quantum chemical calculation is conducted
- `num_relaxation` : Number of optimization during relaxation
- `scale` : **Help Me... lkh**
- `calculator` : quantum chemical calculation program name (currently, `gaussian` and `orca` are available). If you want to use another quantum chemical calculation program, you can implement `new_qc_calculator.py` by referring to `./Calculator/template.py`.
- `unit` : Energy unit (e.g. Hartree, kcal/mol, kJ/mol ...)

**Example Code for diels-alder reaction TS search**
```
export PYTHONPATH=<git clone path>
python run.py \
	-id example/diels \
	-sd example/diels
```

## Output
Output files will be saved in `save directory`. Default option is same with `input directory`

### output.log
Brief description and output logs of `pyMCD`
```
##### Scanning information ######
working_directory: /home/junhkim/ts/pyMCD/examples/diels
command: g16
Energy: Hartree

###### qc_input ######
#N pm6 scf(xqc)

Num relaxation: 5

###### Reactant information ######
charge: 0    multiplicity: 1
16
C 0.0 0.0 0.0
C 0.031015 0.003219 1.336957
C 1.258681 0.002903 2.138678
C 2.495347 -0.002057 1.629726
C 2.266065 -2.426401 0.45058
C 1.151923 -2.426222 -0.276626
H 0.179261 -2.671188 0.128753
H 1.129968 -2.207774 -1.33507
H 2.286373 -2.671063 1.504215
H 3.243833 -2.207765 0.044802
H 3.383219 -0.043833 2.243413
H 2.710598 0.074807 0.572655
H 1.1002 -0.033736 3.219848
H -0.895194 -0.034719 1.916646
H -0.918727 -0.041694 -0.566462
H 0.881426 0.078574 -0.621778

########## Scanning coordinates ##########
(1, 6): 2.7 -> 1.54, -0.1933 angstrom per step
(4, 5): 2.7056 -> 1.54, -0.1943 angstrom per step
Computing node: master

Starting time: 2022-07-25 14:02:04.747527
Set up done!!!
Initialization, Relaxing ....
Relaxation finished! Start scanning ....
[2022-07-25 14:02:08.065] Progress:  (1, 6): 0/6 (4, 5): 1/6
[2022-07-25 14:02:18.950] 0.3920E-2 Hartree has Increased after 6 force calls! 0:00:12.557 Taken ...
...
[2022-07-25 14:04:06.041] Progress:  (1, 6): 6/6 (4, 5): 6/6
[2022-07-25 14:04:11.715] 0.1396E-1 Hartree has Decreased after 6 force calls! 0:00:07.276 Taken ...
[2022-07-25 14:04:11.717077] Scan completed ...
Total 85 force calls performed ...
End time: 2022-07-25 14:04:12.941507
Taken time: 0:02:08.193980

```
### pathway.xyz
Pathway searched by `pyMCD`
```
0
0.0776422942853723
C 0.0 0.0 0.0
C 0.031015 0.003219 1.336957
C 1.258681 0.002903 2.138678
C 2.495347 -0.002057 1.629726
C 2.266065 -2.426401 0.45058
C 1.151923 -2.426222 -0.276626
H 0.179261 -2.671188 0.128753
H 1.129968 -2.207774 -1.33507
H 2.286373 -2.671063 1.504215
H 3.243833 -2.207765 0.044802
H 3.383219 -0.043833 2.243413
H 2.710598 0.074807 0.572655
H 1.1002 -0.033736 3.219848
H -0.895194 -0.034719 1.916646
H -0.918727 -0.041694 -0.566462
H 0.881426 0.078574 -0.621778

...

12
-0.0018082270031678739
C 0.100391 -0.565247 -0.033205
C 0.105619 -0.04404 1.373629
C 1.241504 -0.022778 2.07956
C 2.505306 -0.550895 1.466422
C 2.234325 -1.855663 0.696674
C 0.995067 -1.805221 -0.216555
H 0.386807 -2.712615 -0.034397
H 1.312382 -1.861647 -1.274666
H 2.104802 -2.675434 1.430974
H 3.126446 -2.119199 0.09899
H 3.281606 -0.72351 2.23438
H 2.919783 0.21348 0.775945
H 1.306285 0.35593 3.093065
H -0.841699 0.315363 1.761944
H -0.929587 -0.801621 -0.361092
H 0.454945 0.24615 -0.705328
```
### profile.log
Energy for each step in `pyMCD` pathway
```
Original Energy (Hartree)    Relative Energy (Hartree)
0.0776422942853723   0.0
0.08156201812803866      0.003919723842666367
0.08557842341450855      0.007936129129136257
0.09253669576114133      0.014894401475769037
0.09914460102721953      0.021502306741847235
0.1082233734842189   0.030581079198846603
0.11285269893329185      0.03521040464791955
0.10266648033819747      0.02502418605282518
0.07821241571820944      0.0005701214328371479
0.05345626534766389      -0.02418602893770841
0.026040947061088282     -0.051601347224284014
0.012155719733485644     -0.06548657455188665
-0.0018082270031678739   -0.07945052128854017
Maxima point: 6
```
### profile.png
Relative energy plot

![ex_screenshot](./img/screenshot.png)
### ts.xyz
Geometry of maximal point, `xyz format` with energy
```
16
0.11285269893329185
C 0.111814 -0.215683 -0.038194
C 0.071411 -0.054418 1.327584
C 1.258754 -0.056387 2.096921
C 2.486809 -0.219374 1.496698
C 2.19853 -2.123113 0.602829
C 1.051723 -2.112043 -0.159808
H 0.135785 -2.566401 0.199233
H 1.09879 -1.998238 -1.236483
H 2.209565 -2.583751 1.583162
H 3.175438 -2.03345 0.143827
H 3.379905 -0.378092 2.086275
H 2.691119 0.152435 0.497421
H 1.168918 -0.09791 3.182018
H -0.883297 -0.091713 1.853756
H -0.791863 -0.375329 -0.610443
H 0.940588 0.146519 -0.638458
```
## Citation
If this code was helpful to you, please cite as below
```
Thank you !
```

## References?

## License
