import os
import pytest
import psutil
import itertools
from subprocess import Popen, PIPE, STDOUT

TIMEOUT = 60 * 15  # 15 minute timeout holder

PERCENT = 2  # Random Run 2 % test cases from all combinations

# VCF split testing script

inputFile = [['--vcf', 'merged_chr1_10000.vcf.gz']]
outPrefix = [['--out-prefix', 'merge_chr1_10000']]
phaseAlgo = [['--phase-algorithm beagle'], ['--phase-algorithm shapeit']]

overwrite = [['--overwrite']]


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
        ids.append('python vcf_split.py ' + ' '.join(i))
    return ids


filterArgs0 = list(itertools.product(inputFile, modelFile, model, outPrefix, splitFile, splitMethod, outFormat, overwrite))
filterArgs1 = list(itertools.product(inputFile, model, outPrefix, splitFile, splitMethod, outFormat, Remaining_options, overwrite))

all_combinations = list_of_lists(filterArgs0) + list_of_lists(filterArgs1)


@pytest.mark.parametrize('split_args', all_combinations, ids=id_func(all_combinations))
@pytest.mark.timeout(TIMEOUT)      # 15 minutes test case timeout
def test_(split_args):
    p = Popen(['python ' + os.getcwd() + '/vcf_split.py ' + ' '.join(flatten(split_args))], shell=True, stdin=PIPE,
              stdout=PIPE, stderr=STDOUT, close_fds=True)

    q = psutil.Process(p.pid)
    mem = q.memory_info()[0] / 2. ** 10
    cpu = q.cpu_percent(interval=0.1)
    print('Memory Usage:', mem, 'kB')
    print('CPU Usage:', cpu, '%')

    output = p.stdout.read().decode().split('\n')
    p.communicate()
    assert p.returncode == 0, output


# pytest --html=reports/vcf_split_test_suite.html --capture=fd --self-contained-html --capture=sys --tb=no test_4.py
