from setuptools import setup, find_packages

setup(name="scPathoQuant",
      version="1.2.0",
      description="The goal of this package is to accurately align and quantify viral reads for 10x single cell data. \
                   This software integrates viral counts and viral gene counts into 10x files (features.tsv.gz and matrix.mtx.gz \
                  in the filtered_feature_bc_matrix folder) so that softwares such as seurat can be used to analyze data",
      author=["Leanne Whitmore"],
      author_email=["leanne382@gmail.com"],
      platforms=["linux"],
      keywords='single-cell, viral/pathogen, quantification',
      packages=find_packages(),
      install_requires=[
        "pandas==2.0.3",
        "argparse==1.4.0",
        "htseq==2.0.5",
        "scipy==1.10.1",
        "seaborn==0.13.0",
        "pysam==0.22.0",
        "setuptools>=59.0.0",
        "matplotlib==3.3.2",
        "GFFUtils==0.12",
        "pyGenomeTracks==3.6"
      ],
      include_package_data=True,
      zip_safe=False,
      scripts=[
          'scpathoquant',
      ]
)