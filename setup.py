import os
import re
import sys
from glob import glob
import shutil
from setuptools import setup, find_packages, Extension
from setuptools import Distribution
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
import subprocess
from pathlib import Path

__pkg_name__ = 'pepper'
__author__ = 'Kishwar Shafin'
__description__ = 'RNN based assembly polisher.'


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
            if cmake_version < '3.11':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")
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
            self.move_output(ext)

    def move_output(self, ext):
        build_temp = Path(self.build_temp).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve().parents[0] / "build" / self.get_ext_filename(ext.name)
        source_path = build_temp / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parents[0]
        dest_directory.mkdir(parents=True, exist_ok=True)
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


if __name__ == '__main__':
    # check python3 version
    pymajor, pyminor = sys.version_info[0:2]
    if (pymajor < 3) or (pymajor <= 3 and pyminor < 5):
        raise RuntimeError(
            'PEPPER requires 3.5 higher.')
    python_dependencies = get_dependencies()
    setup(
        name='PEPPER',
        version='1.0',
        packages=['', 'modules/python', 'modules/python/models', 'modules/python/helper'],
        package_dir={'modules/python': 'modules/python',
                     'modules/python/models': 'modules/python/models',
                     'modules/python/helper': 'modules/python/helper'},
        url='https://github.com/kishwarshafin/pepper',
        author=__author__,
        description=__description__,
        python_requires='>=3.5.*',
        install_requires=python_dependencies,
        entry_points={
            'console_scripts': [
                '{0} = pepper:main'.format(__pkg_name__),
                '{0}_train = pepper_train:main'.format(__pkg_name__),
            ]
        },
        ext_modules=[CMakeExtension('PEPPER')],
        cmdclass={
            'build_ext': CMakeBuild
        },
        zip_safe=False,
    )
