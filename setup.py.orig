import os
import re
import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
import subprocess

__pkg_name__ = 'pepper'
__author__ = 'Kishwar Shafin'
__description__ = 'RNN based genome inference tool.'


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                                   out.decode()).group(1))
            print("CMAKE VERSION: ", cmake_version)
            if cmake_version < '3':
                raise RuntimeError("CMake >= 3 is required.")
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        build_directory = os.path.abspath(self.build_temp)

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + build_directory,
            '-DPYTHON_EXECUTABLE=' + sys.executable
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        # Assuming Makefiles
        build_args += ['--', '-j2']

        self.build_args = build_args

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # CMakeLists.txt is in the same directory as this setup.py file
        cmake_list_dir = os.path.abspath(os.path.dirname(__file__))
        print('-'*10, 'Running CMake prepare', '-'*40)
        subprocess.check_call(['cmake', cmake_list_dir] + cmake_args,
                              cwd=self.build_temp, env=env)

        print('-'*10, 'Building extensions', '-'*40)
        cmake_cmd = ['cmake', '--build', '.'] + self.build_args
        subprocess.check_call(cmake_cmd,
                              cwd=self.build_temp)

        # Move from build temp to final position
        for ext in self.extensions:
            print(ext)
            self.move_output(ext)

    def move_output(self, ext):
        source_path = os.path.abspath(self.build_temp + "/" + self.get_ext_filename(ext.name))
        print(source_path)
        dest_directory = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name))) + "/" + str(ext.name).lower() + "/build/"
        print(dest_directory)
        os.makedirs(dest_directory, exist_ok=True)

        dest_path = dest_directory + self.get_ext_filename(ext.name)

        self.copy_file(source_path, dest_path)


def get_dependencies():
    # create requirements from requirements.txt
    dir_path = os.path.dirname(__file__)
    install_requires = []
    with open(os.path.join(dir_path, 'requirements.txt')) as fh:
        reqs = (
            r.split('#')[0].strip()
            for r in fh.read().splitlines() if not r.strip().startswith('#')
        )
        for req in reqs:
            if req.startswith('git+https'):
                req = req.split('/')[-1].split('@')[0]
            install_requires.append(req)
    return install_requires


def get_version():
    version = {}
    with open("./pepper/version.py") as fp:
        exec(fp.read(), version)
    return version['__version__']


def get_long_description():
    from os import path
    this_directory = path.abspath(path.dirname(__file__))
    kwargs = {'encoding':'utf-8'} if sys.version_info.major == 3 else {}
    with open(path.join(this_directory, 'README.md'), **kwargs) as f:
        long_description = f.read()
    long_description_content_type = 'text/markdown'

    return long_description, long_description_content_type


if __name__ == '__main__':
    # check python3 version
    pymajor, pyminor = sys.version_info[0:2]
    if (pymajor < 3) or (pymajor <= 3 and pyminor < 5):
        raise RuntimeError(
            'PEPPER requires 3.5 higher.')
    python_dependencies = get_dependencies()
    __long_description__, __long_description_content_type__ = get_long_description()

    setup(
        name='pepper_polish',
        version=get_version(),
        packages=['pepper/', 'pepper/modules/python', 'pepper/modules/python/models', 'pepper/modules/python/helper',
                  'pepper_variant/', 'pepper_variant/modules/python', 'pepper_variant/modules/argparse', 'pepper_variant/modules/python/models'],
        package_dir={'pepper': 'pepper',
                     'pepper/modules/python': 'pepper/modules/python',
                     'pepper/modules/python/models': 'pepper/modules/python/models',
                     'pepper/modules/python/helper': 'pepper/modules/python/helper',
                     'pepper_variant': 'pepper_variant',
                     'pepper_variant/modules/python': 'pepper_variant/modules/python',
                     'pepper_variant/modules/argparse': 'pepper_variant/modules/argparse',
                     'pepper_variant/modules/python/models': 'pepper_variant/modules/python/models',
                     'pepper_variant/modules/python/helper': 'pepper_variant/modules/python/helper'
                     },
        url='https://github.com/kishwarshafin/pepper',
        author=__author__,
        author_email="kishwar.shafin@gmail.com",
        description=__description__,
        long_description="Please visit (GitHub)[https://github.com/kishwarshafin/pepper] for description",
        long_description_content_type=__long_description_content_type__,
        python_requires='>=3.5.*',
        install_requires=python_dependencies,
        entry_points={
            'console_scripts': [
                '{0} = {0}.{0}:main'.format(__pkg_name__),
                '{0}_train = {0}.{0}_train:main'.format(__pkg_name__),
                '{0}_variant = {0}_variant.{0}_variant:main'.format(__pkg_name__),
                '{0}_variant_train = {0}_variant.{0}_variant_train:main'.format(__pkg_name__),
            ]
        },
        ext_modules=[CMakeExtension('PEPPER'), CMakeExtension('PEPPER_VARIANT')],
        cmdclass={
            'build_ext': CMakeBuild
        },
        zip_safe=False,
    )
