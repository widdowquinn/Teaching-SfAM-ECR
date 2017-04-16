FROM andrewosh/binder-base

USER root

# Add dependency
RUN apt-get update

# Install BLAST+
RUN apt-get install -y ncbi-blast+

USER main

# Get rid of token/password requests
RUN mkdir -p $HOME/.jupyter
RUN echo "c.NotebookApp.token = ''" >> $HOME/.jupyter/jupyter_notebook_config.py

# Install requirements for Python 3
RUN /home/main/anaconda/envs/python3/bin/pip install -r requirements.txt

# Install ipywidgets under conda
RUN conda config --add channels conda-forge
RUN conda install -y ipywidgets widgetsnbextension

# Enable widgets extension
RUN /home/main/anaconda/envs/python3/bin/jupyter nbextension enable --py --sys-prefix widgetsnbextension

# Add new kernel
RUN /home/main/anaconda/envs/python3/bin/python -m ipykernel install --user --name Python3_ibioic_course --display-name "Python 3 (SfAM)"
