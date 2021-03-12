import admixture
import os
import pytest
import psutil
import itertools
from subprocess import Popen, PIPE, STDOUT

# Admixture function testing script
inputFile = [['--file', 'hapmap3.bed']]
population = [['--pop', '4']]
randSeed = [None, ['--random-seed', '12345']]
method = [None, ['--method', 'em']]
acceleration = [None, ['--acceleration', 'sqs5'], ['--acceleration', 'qn5']]
convergence = [None, ['--major-converge', '2'], ['--minor-converge', '2']]
bootstrap = [None, ['--bootstrap', '200']]


def flatten(tuple_of_lists):
    nested_list = list(tuple_of_lists)
    rt = []
    for i in nested_list:
        if isinstance(i, list):
            rt.extend(flatten(i))
        else:
            rt.append(i)
    return list(filter(None.__ne__, rt))


def list_of_lists(list_of_tuple):
    flat_list = []
    for i in list_of_tuple:
        flat_list.append(flatten(i))
    return flat_list


def id_func(args):
    ids = []
    for i in args:
        ids.append('python admixture.py ' + ' '.join(i))
    return ids

filterArgs = list(itertools.product(inputFile, population, randSeed, method, acceleration, bootstrap))
# filterArgs1 = list(itertools.product(inputFileList, outFormatList, outPrefixList, list_of_lists(itertools.combinations(filter(None.__ne__, Remaining_options1), 2)), overwrite))

all_combinations = list_of_lists(filterArgs)


@pytest.mark.parametrize('admixture_args', all_combinations, ids=id_func(all_combinations))
def test_admixture(admixture_args):
    p = Popen(['python ' + os.getcwd() + '/admixture.py ' + ' '.join(flatten(admixture_args))], shell=True, stdin=PIPE,
              stdout=PIPE, stderr=STDOUT, close_fds=True)

    q = psutil.Process(p.pid)
    mem = q.memory_info()[0] / 2. ** 10
    cpu = q.cpu_percent(interval=0.1)
    print('Memory Usage:', mem, 'kB')
    print('CPU Usage:', cpu, '%')

    output = p.stdout.read().decode().split('\n')
    p.communicate()
    assert p.returncode == 0, output


# pytest --html=reports/admixture_test_suite.html --capture=fd --self-contained-html --capture=sys --tb=no test_6.py