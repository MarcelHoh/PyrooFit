PyrooFit
========

PyrooFit is a fit framework for python and pandas DataFrames on top of the ROOT.RooFit package developed by @simonUU.

The package allows for simple fits of standard PDFs and easy setup of custom PDFs in one or more fit dimensions. This fork replaces the default root style plotting with python matplotlib based plotting.

Example
-------

Simple fit and plot of a Gaussian Distribution:

```python
from pyroofit.models import Gauss
import numpy as np

data = np.random.normal(0, 1, 1000)

pdf = Gauss(('x', -3, 3), mean=(-1, 0, 1))
pdf.fit(data)
pdf.plot('example_gauss.pdf',)

pdf.get()

```

A more complex example on combination of Gauss pdf for signal and Polynomial for background:

```python
from pyroofit.models import Gauss, Chebychev
import numpy as np
import pandas as pd
import ROOT



df = {'mass': np.append(np.random.random_sample(1000)*10 + 745, np.random.normal(750, 1, 1000))}
df = pd.DataFrame(df)

x = ROOT.RooRealVar('mass', 'M', 750, 745, 755, 'GeV')  # or x = ('mass', 745, 755)

pdf_sig = Gauss(x, mean=(745, 755), sigma=(0.1, 1, 2), title="Signal")
pdf_bkg = Chebychev(x, n=1, title="Background")

pdf = pdf_sig + pdf_bkg

pdf.fit(df)
pdf.plot('example_sig_bkg.pdf', legend=True)
pdf.get()

```

<img src="http://desy.de/~hohmann/example_sig_bkg.png" width="400" height="400">


Observables can be initialised by a list or tuple with the column name / variable name as first argument, followed
by the range and/or with the initial value and range:
```
x = ('x', -3, 3)
x = ('mass', -3, 0.02, 3)
```

Parameters are initialised with a tuple: `sigma=(0,1)`
or again including a starting parameter: `sigma=(0.01, 0, 1)`
The order here is not important.

All parameters and observables can also be initialised by a `ROOT.RooRealVar`.

Styling
============

The styling is deliberately left to a minimum to allow users to use their own style files. It is recommended to use a stylesheet like those provided by matplotlib or within [b2plot](https://github.com/simonUU/b2plot) (also by @simunUU).

In addition to the styling performed by the stylesheet or modifying `rcParams` there are simple options for `fill`, `hatch` and `color` attached to each pdf and each element of the plot can be individually edited by supplying a dictionary of keyword arguments. For example:

```python
from pyroofit.models import Gauss, Chebychev
import numpy as np
import pandas as pd
import ROOT

plt.style.use('belle2_wg1')

df = {'mass': np.append(np.random.random_sample(1000)*10 + 745, np.random.normal(750, 1, 1000))}
df = pd.DataFrame(df)

x = ROOT.RooRealVar('mass', 'M', 750, 745, 755, 'GeV')  # or x = ('mass', 745, 755)

pdf_sig = Gauss(x, mean=(745, 755), sigma=(0.1, 1, 2), title="Signal",
                color='xkcd:coral', fill=True, fill_alpha=0.3)
pdf_bkg = Chebychev(x, n=1, title="Background",
                color='xkcd:sea blue', hatch='\\')

pdf = pdf_sig + pdf_bkg

pdf.fit(df)
pdf.plot('example_sig_bkg_2.pdf',
         title = 'Style: belle2_wg1',
         title_kwargs = {'fontsize' : 20},
         legend_kwargs = {'markerfirst' : False})
pdf.get()

```
<img src="http://desy.de/~hohmann/example_sig_bkg_2.png" width="400" height="400" alt="use this one">
<img src="http://desy.de/~hohmann/example_sig_bkg_3.png" width="400" height="400">
<img src="http://desy.de/~hohmann/example_sig_bkg_4.png" width="400" height="400">


Installation
============

Dependencies: ROOT (with PyRoot enabled)


* Download this repository

* (recommended) Use or install anaconda python environment

* Activate ROOT installation with python support

* run ``python setup.py install`` in this folder

* run ``python setup.py docs`` to create the documentation

If you do not have your own python installation you can use:
```
python setup.py install --user
PATH=$PATH~/.local/bin
```
If there are still missing packages you might need to install them via
`pip install package --user`.



Planned Features
================

- Improve documentation
- Save and load PDF as yaml
- Plotting in matpltotlib
