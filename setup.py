from setuptools import Extension, setup

module = Extension("spkmeansmodule",
                  sources=[
                    'spkmeans.c',
                    'spkmeansmodule.c',
                    'spkmeans.c',
                    'wam.c',
                    'ddg.c',
                    'lnorm.c',
                    'jacobi.c',
                    'utils.c',
                    'normalizedKEigenvectorsMatrix.c'
                  ])
setup(name='spkmeansmodule',
     version='1.0',
     description='Normalized Spectral Clustering',
     ext_modules=[module])
