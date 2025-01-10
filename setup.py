if __name__ == "__main__":
    from setuptools import setup

    setup(
        name="genomediff",
        version="0.4.0",
        packages=["genomediff"],
        url="https://github.com/Hocnonsense/genomediff-python",
        license="MIT",
        author="Aoran Hu",
        author_email="hwrn.aou@sjtu.edu.cn",
        description="GenomeDiff (*.gd) file reader",
        long_description="GenomeDiff file reader",
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Environment :: Other Environment",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )
