[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pmitev/jupyter/master?filepath=IJQC.116.db%2FIJQC.116.properties.ipynb)
# Folder contents:

 - **[IJQC.116.db](https://github.com/pmitev/jupyter/tree/master/IJQC.116.db)** - contains an example Jupyter notebook to reproduce the results published in "*H-bond and electric field correlations for water in highly hydrated crystals*", Anik Sen, Pavlin D. Mitev, Anders Eriksson, and Kersti Hermansson, International Journal of Quantum Chemistry, **116**, ( 2016), 67-80 **DOI**: [10.1002/qua.25022](http://onlinelibrary.wiley.com/doi/10.1002/qua.25022/abstract)

---

[ASE](https://wiki.fysik.dtu.dk/ase/index.html) [format](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) database with all structures, Bader and Wannier analysis, and OH stretching frequency published in the article.

[Download](https://github.com/pmitev/jupyter/raw/master/IJQC.116.db/IJQC.116.db)

    $> ase db IJQC.116.db -c id,name,method,pbc,charge,volume,energy
    id|name         |method|pbc|charge|  volume|  energy
     1|Al(NO3)3.9H2O|PBE-D2|TTT| 0.000|1421.939|-882.352
     2|MgSO4.11H2O  |PBE-D2|TTT| 0.000| 680.318|-405.280
     3|MgSO4.7H2O   |PBE-D2|TTT| 0.000| 933.902|-571.405
     4|Na2CO3.10H2O |PBE-D2|TTT| 0.000|1219.089|-752.414

You can use web browser to look on common prperties

    $> ase db IJQC.116.db -w
     * Running on http://0.0.0.0:5000/ (Press CTRL+C to quit)

![ ](http://www.teoroo.kemi.uu.se/wp-content/uploads/2017/11/ase-s.png)


Command line interface:

    $> ./ASE_db_tools.py -h
    usage: ASE_db_tools.py [-h] [-id]
                           (--add-VASP-tructure | --add-VASP-bader | --add-VASP-wannier | --add-VASP-DDEC6 | --add-symmetry | --add-OH-PES | --add-DVR-results | --prop-H2O-bader-dipoles | --prop-H2O-wannier-dipoles | --prop-H2O-DDEC6-dipoles | --EF-H2O-bader | --EF-H2O-wannier | --EF-H2O-DDEC6)
                           [-d] [-f] [-l] [-i] [-s] [-chg] [-v]
                           dbname
    
    ASE database tools
    
    positional arguments:
      dbname                database name
    
    optional arguments:
      -h, --help            show this help message and exit
      -id                   number of the row in the database
      --add-VASP-tructure   Add VASP structre to the database
      --add-VASP-bader      Add VASP Bader calculation results to the database
      --add-VASP-wannier    Add VASP Wannier calculation results to the database
      --add-VASP-DDEC6      Add VASP DDEC6 calculation results to the database
      --add-symmetry        Find and add symmetry for the structure in the
                            database
      --add-OH-PES          Add OH stretch PES in the database
      --add-DVR-results     Add DVR results for the sytructure in the database
      --prop-H2O-bader-dipoles
                            Calculate H2O Bader dipoles for the structure in the
                            database
      --prop-H2O-wannier-dipoles
                            Calculate H2O Wannier dipoles for the structure in the
                            database
      --prop-H2O-DDEC6-dipoles
                            Calculate H2O DDEC6 dipoles for the structure in the
                            database
      --EF-H2O-bader        Calculate Bader EF at H positions for H2O molecules
                            for the structure in the database
      --EF-H2O-wannier      Calculate Wannier EF at H positions for H2O molecules
                            for the structure in the database
      --EF-H2O-DDEC6        Calculate DDEC6 EF at H positions for H2O molecules
                            for the structure in the database
      -d , --directory      Specify folder input
      -f , --file           Specify file input
      -l , --label          Specify label input
      -i , --index          Specify index input
      -s , --symprec        Specify symmetry precision input
      -chg , --core_charges 
                            Specify core_charges in json format. '{ "H" : 1.0 ,
                            "O" : 6.0 }'
      -v, --verbose         Verbose output 

Example:

    $ ./ASE_db_tools.py  IJQC.116.db -id 1 --prop-H2O-wannier-dipoles --core_charges '{"H" : 1, "O": 6 , "W": -2 , "Al": 3 , "N": 5}' | head
      Origin at atom idx: 88  core_charge: 6
    ============================================
    idx:  213   Chg:   -2   dist: 0.478127
    idx:  286   Chg:   -2   dist: 0.293650
    idx:    5   Chg:    1   dist: 1.000710
    idx:  302   Chg:   -2   dist: 0.479148
    idx:  337   Chg:   -2   dist: 0.359560
    idx:    4   Chg:    1   dist: 0.999393
    ============================================
      Total Dipole= 3.9117 [Debye]
