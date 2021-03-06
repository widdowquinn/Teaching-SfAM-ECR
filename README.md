# SfAM ECS 2017 Presentation and Workshop

## Table of Contents

1. [Introduction](#introduction)
2. [Schedule](#schedule)
3. [Presentations](#presentations)
4. [Workshop Notebooks](#notebooks)
5. [Run Workshop in the Cloud on MyBinder](#mybinder)
6. [Run Workshop on Your Own Machine](#local)


<a id="introduction"></a>
## Introduction
Welcome to the "Not All Bioinformatics is NGS" part of the SfAM ECS symposium. This page describes the contents of the repository, provides useful links, and instructions for how to reproduce and run the workshop materials on your own machine, or run them in the cloud.

<a id="schedule"></a>
## Schedule
* [[poster]](posters/SfAM_Bioinformatics_A4.pdf)

| Time          | Activity                                         |
| :-----------: | :----------------------------------------------- |
| 10:00 - 11:00 | Registration                                     |
| 10:30 - 11:00 | Posters & Refreshments                           |
| 11:00 - 12:00 | Student oral presentations (4x 15 min slots)     |
| 12:00 - 13:00 | Introduction to Bioinformatics (2x 30 min slots) |
| 13:00 - 14:00 | Posters, Exhibitions, and Lunch                  |
| 14:00 - 15:30 | Workshop 1                                       |
| 15:30 - 16:00 | Posters & Refreshments                           |
| 16:00 - 17:30 | Workshop 2                                       |
| 17:30 - 18:00 | Conference Close (Prizes & Networking)           |

<a id="presentations"></a>
# Presentations

The links below will take you to the presentations in an interactive form, in your browser.

* [My Life in Sequences](https://widdowquinn.github.io/Teaching-SfAM-ECS/presentation/slides00-keynote.html)
* [Reproducible Research](https://widdowquinn.github.io/Teaching-SfAM-ECS/presentation/slides01-reproducible_research.html)


<a id="notebooks"></a>
# Workshop Notebooks

The links below will take you to the workshop notebooks in static form at `GitHub`, in your browser.

* [**00 - Introduction to Jupyter**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/00-jupyter_intro.ipynb) 
* [**01 - Thinking Statistically**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/01-thinking_statistically.ipynb)
  * [**01a - Correlations**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/01a-correlations.ipynb)
  * [**01b - Classifiers**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/01b-clasifiers.ipynb)
* [**02 - Automation**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/02-automation.ipynb)
  * [**02a - Python Basics**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/02a-python.ipynb)
  * [**02b - Automating BLAST**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/02b-blast.ipynb)
* [**03 - Integrating Data**](https://github.com/widdowquinn/Teaching-SfAM-ECS/blob/master/workshop/03-integrating_data.ipynb)

<a id="mybinder"></a>
## Run Workshop in the Cloud on MyBinder

MyBinder is a service that allows you to run Jupyter notebooks in the cloud (i.e. on someone else's hardware). This course is provided as interactive notebooks *via* this service, and you can start up an instance by clicking on the button below:

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/widdowquinn/teaching-sfam-ecs)

This approach has the advantage that it requires no further configuration: all materials and tools are pre-installed, and ready to go.


<a id="local"></a>
## Run Workshop on Your Own Machine

Setting up the course materials to run on your own machine involves following stages, which will be described one-by-one below:

1. Clone the repository
2. Ensure the bioinformatics package prerequisites are installed
3. Create a virtual environment for the Jupyter notebooks

To install and run these materials locally, you will require:

* The Anaconda Python distribution [[download]](https://www.continuum.io/downloads)

### Download the materials

The course materials (notebooks, presentations, data) are available for download at the link below:

* [https://github.com/widdowquinn/Teaching-SfAM-ECS/archive/v1.0.zip](https://github.com/widdowquinn/Teaching-SfAM-ECS/archive/v1.0.zip)

Download the compressed (`.zip`) file to a suitable location and unpack it. This will place the materials in a subdirectory called `Teaching-SfAM-ECR`. 

Alternatively, if you have [`git`](https://git-scm.com/) installed, you can clone the most current version of the materials to your local machine with the command:

```bash
git clone https://github.com/widdowquinn/Teaching-SfAM-ECS.git
```

### Run the materials

* Open a terminal on your machine
* Navigate to a convenient location and create a new directory for the course materials (e.g. `sfam-workshop`)

```bash
mkdir sfam-workshop
cd sfam-workshop
```

* Download or clone the course materials to this location (see above) and extract them, if necessary
* Change directory to the course materials

```bash
cd Teaching-SfAM-ECS        # if you cloned the repository
cd Teaching-SfAM-ECS-1.0    # if you extracted the .zip/.tar.gz file
```

* Add a new `Anaconda` channel for `bioconda`

```bash
conda config --add channels bioconda
```

* Install `Jupyter` in the `conda` environment

```bash
conda install jupyter
```

* Create a new `conda` environment called `sfam`, installing the necessary `Python` libraries (accept the changes requested).

```bash
conda create --name sfam python=3 jupyter numpy pandas seaborn bioservices biopython
```

* Activate the `conda` environment:

On Windows:

```bash
activate sfam
```

On OSX/Linux

```bash
source activate sfam
```

You should see that your command-line prompt has changed, and now starts with `(sfam)`. This indicates that you are now working within a 'protected' `Python` environment, specific for this workshop.

* Add a new kernel to `Jupyter`, and enable iPython widgets

```bash
python -m ipykernel install --user --name Python3_SfAM --display-name "Python 3 (SfAM)"
jupyter nbextension enable --py --sys-prefix widgetsnbextension
```

* Start the notebooks with `jupyter notebook index.ipynb` to begin at the top page of the workshop

```bash
jupyter notebook index.ipynb
```



