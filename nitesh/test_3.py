import os
import pytest
import psutil
import itertools
from random import shuffle
from subprocess import Popen, PIPE, STDOUT

TIMEOUT = 60 * 15  # 15 minute timeout holder

PERCENT = 2  # Random Run 2 % test cases from all combinations

# VCF stat sampler testing script

statFile = [['--statistic-file merged_chr1_10000.windowed.weir.fst'], ['--statistic-file merged_chr1_10000.Tajima.D'], ['--statistic-file merged_chr1_10000.windowed.pi']]
outPrefix = [['--out-prefix', 'merged_chr1_10000']]
calcStat = [['--calc-statistic windowed-weir-fst'], ['--calc-statistic TajimaD'], ['--calc-statistic window-pi']]
sampling = [None, ['--sampling-scheme uniform'], ['--sampling-scheme random']]
uniformBins = [None, ['--uniform-bins', '5']]
sampleSize = [None, ['--sample-size', '20']]
randomSeed = [None, ['--random-seed', '100']]
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
        ids.append('python stat_sampler.py ' + ' '.join(i))
    return ids


def random_tests_select(list_a, percentage):
    shuffle(list_a)
    percentage = percentage / 100
    count = int(len(list_a) * percentage)
    if not count: return []  # edge case, no elements removed
    list_a[-count:], list_b = [], list_a[-count:]
    return list_b


calcArgs0 = list(itertools.product(statFile[0], calcStat[0], outPrefix, sampling, uniformBins, sampleSize, randomSeed, overwrite))
calcArgs1 = list(itertools.product(statFile[1], calcStat[1], outPrefix, sampling, uniformBins, sampleSize, randomSeed, overwrite))
calcArgs2 = list(itertools.product(statFile[2], calcStat[2], outPrefix, sampling, uniformBins, sampleSize, randomSeed, overwrite))

all_combinations = list_of_lists(calcArgs0) + list_of_lists(calcArgs1) + list_of_lists(calcArgs2)


@pytest.mark.parametrize('stat_sampler', all_combinations, ids=id_func(all_combinations))
@pytest.mark.timeout(TIMEOUT)      # 15 minutes test case timeout
def test_(stat_sampler):
    p = Popen(['python ' + os.getcwd() + '/stat_sampler.py ' + ' '.join(flatten(stat_sampler))], shell=True, stdin=PIPE,
              stdout=PIPE, stderr=STDOUT, close_fds=True)

    q = psutil.Process(p.pid)
    mem = q.memory_info()[0] / 2. ** 10
    cpu = q.cpu_percent(interval=0.1)
    print('Memory Usage:', mem, 'kB')
    print('CPU Usage:', cpu, '%')

    output = p.stdout.read().decode().split('\n')
    p.communicate()
    assert p.returncode == 0, output


# pytest --html=reports/vcf_stat_sampler_test_suite.html --capture=fd --self-contained-html --capture=sys --tb=no test_3.py
