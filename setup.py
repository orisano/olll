from setuptools import setup


setup(
    name="olll",
    version="1.0.2",
    author="Nao Yonashiro",
    author_email="owan.orisano@gmail.com",
    description="A Python3 implementation of the LLL",
    license="MIT",
    url="https://github.com/orisano/olll",
    py_modules=["olll"],
    install_requires=["numpy"],
    extras_require={
        "develop": ["pytest"]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
)
