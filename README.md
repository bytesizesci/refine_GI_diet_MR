# Dietary Intake Mendelian Randomization

This project corresponds to our manuscript, *Dietary Intake Mendelian Randomization: Assessment and Development of Methods for Instrument Selection and Robust Inference*, written by 
Kristen J. Sutton, Julie Gervis, Moomal Jatoi, Liang-Dar Hwang, Audrey Hendricks, Debashis Ghosh, Kenneth Westerman, and Joanne B. Cole.

Read our preprint: [https://doi.org/10.1101/2025.06.26.25330002](https://doi.org/10.1101/2025.06.26.25330002 "Awesome dietary intake Mendelian randomization paper") 

## 🗂️ Project Structure

```text
project-root/
├── scripts_archive/              # Archive of analysis and pipeline scripts
├── results_archive/              # Archive of figures and tables
├── working_example/              # Under development
│   ├── data/                     # Contains data to run the working_example
│   │   ├── reference_1KGP3_HG19/ # Contains reference data
│   └── scripts/                  # Contains the scripts to run the working_example
├── .gitignore
└── README.md
```

## :eyeglasses: Prerequisites 

* Download the [Neale Lab's GWAS round 2](https://www.nealelab.is/uk-biobank) summary statistic files. These files are used for the PheWAS multivariable MR (MVMR) pipeline.
* Run [PheWAS-based clustering](https://www.nature.com/articles/s41467-024-45655-8) on your exposure and outcome of interest to get a list of suspicious SNPs and traits used in the PheWAS MVMR pipeline. PheWAS-based clustering GitHub can be found [here](https://github.com/LizaDarrous/PheWAS-cluster). Note, we call this method PheWAS-ttest in our manuscript.

## 🧬 Data availability

Data used in the manuscript are publicly and freely available without restriction at the GWAS catalog (Supplemental Table 1: LDL-C, GCST90239658; CVD, GCST90132314; TG, GCST90239664; ALT, GCST90013405; Cirrhosis of liver, GCST90013405; Height, GCST006901) and the Type 2 Diabetes Knowledge Portal by the dataset name, UK Biobank dietary habit GWAS.
