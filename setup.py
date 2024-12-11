from setuptools import setup


#jbw 2024, here added four required pacakges!!
#scikit-learn, requests','fpdf','pillow'
setup(name='abc4pwm',
      version='0.1',
      description='Affinity Based Clustering for Position Weight Matrices',
      url='',
      author='Omer Ali',
      author_email='omerali.0191@gmail.com',
      include_package_data=True,
      license='MIT',
      packages=['abc4pwm'],
      install_requires=[
          'numpy',
          'scipy',
          'weblogo',
          'scikit-learn',
          'lxml','requests','fpdf','pillow',
          'bs4'
      ],
    scripts=['abc4pwm/abc4pwm']
      )
