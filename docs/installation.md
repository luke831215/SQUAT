# Installation
---

The installation of SQUAT is very simple.  
If you are getting a permission error during the installation, consider running with `sudo`, or create a virtual python environment.

## Docker installation
hub: https://hub.docker.com/r/luke831215/squat/

	docker pull luke831215/squat
	docker run --name=squat_docker -it luke831215/squat

-- Running container --

	git clone https://github.com/luke831215/SQUAT
	cd SQUAT; make install

## Github installation

### Dependency
- [Python3](https://www.python.org/downloads/)
- [Matplotlib](https://matplotlib.org/faq/installing_faq.html#python-org-python)
- [Numpy](https://askubuntu.com/questions/765494/how-to-install-numpy-for-python3)
- [Beautiful Soup (+ lxml parser)](https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup) 

`pip3 install matplotlib numpy beautifulsoup4 lxml`

First, you need to clone down the repository.

	git clone https://github.com/luke831215/squat

Then, go to the SQUAT directory and execute the makefile.
	
	cd SQUAT; make install

After the installation, start running the tool with **squat.sh**.

## Example data

For trial purposes, we extract 25000 reads from the specie **Saccharomyces cerevisiae** and its assembly. Run the following command to start using. See [Usage](usage.md) and [Output](output.md) for more details.

	./squat.sh example/SEQ.fastq -o example -r example/ASSEMBLY.fasta

After finishing, open `SEQ.html` in example directory to begin.