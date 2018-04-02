# Installation

## Dependency
- Python3
- Matplotlib
- Numpy

## Github

The installation of SQUAT is very simple.  
If you are getting a permission error during the installation, consider running with `sudo`, or create a virtual python environment.

First, you need to clone down the repository.

	git clone https://github.com/luke831215/squat

Then, go to the SQUAT directory and execute the makefile.
	
	cd SQUAT; make install

After the installation, start running the tool with **squat.sh**.

## Example

As example data, we extract 25000 reads from the specie **Saccharomyces cerevisiae** and its assembly. Run the following command to start using. See [Usage](usage.md) and [Output](output.md) for more details.

	bash squat.sh example/SEQ.fastq -o example -r example/ASSEMBLY.fasta

After finishing, open `SEQ.html` in example directory to begin.