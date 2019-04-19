# 'xlrd' Package Lab

This lab provides some practice working with library code in Python: 
* creating a virtual environment (a package-level namespace you can manage seperately from your global set of packages) 
* using the packages [`xlrd`](https://xlrd.readthedocs.io/en/latest/) for reading and writing Microsoft Excel spreadsheets and [`statistics`](https://docs.python.org/3/library/statistics.html) to describe the data you read.

<hr>

#### Step-by-step to get started

* On the Linux lab machines, clone this repo (and remove the `.git` directory if you've cloned it inside your grading repo).

* Inside this repo, create a python3 virtual environment named `venv` using the tool `virtualenv`:    
```virtualenv -p python3 venv ```

* Activate the `venv` virtual environment.
```source venv/bin/activate```

* Use `pip3` to install the packages `xlrd` and `statistics` into this virtual environment:     
```pip3 install xlrd```    
```pip3 install statistics```    

* Follow the instructions inside `main.py` to complete the lab. You're welcome to change the code however you'd like, but you should only need to add code in the sections marked `[FILL IN HERE]`. All of the functions outside of the main block are there to help you with debugging or manipulating the data in order to answer the two questions in the main block.

* Deactivate your virtual environment:    
```deactivate```
