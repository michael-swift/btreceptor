import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="btreceptor", # Replace with your own username
    version="0.0.1",
    author="Derek Croote",
    author_email="author@example.com",
    description="bt receptor package for blue tooth receptor analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dcroote/btreceptor.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data=True
)
