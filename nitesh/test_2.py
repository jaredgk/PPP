import os
import pytest
import psutil
import itertools
from random import shuffle
from subprocess import Popen, PIPE, STDOUT

TIMEOUT = 60 * 15  # 15 minute timeout holder

PERCENT = 2  # Random Run 2 % test cases from all combinations

# VCF calc testing script

inputFileList = [['--vcf', 'merged_chr1_10000.vcf.gz']]
modelFile = [['--model-file', 'input.model']]
model = [['--model', '2Pop']]
outPrefix = [['--out-prefix', 'merged_chr1_10000']]
calcStat = [['--calc-statistic windowed-weir-fst'], ['--calc-statistic weir-fst'], ['--calc-statistic TajimaD'], ['--calc-statistic window-pi'], ['--calc-statistic site-pi'], ['--calc-statistic freq'], ['--calc-statistic het-fit'], ['--calc-statistic het-fis']]
windowSize = [['--statistic-window-size', '10000']]
windowStep = [['--statistic-window-step', '20000']]
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
        ids.append('python vcf_calc.py ' + ' '.join(i))
    return ids


def random_tests_select(list_a, percentage):
    shuffle(list_a)
    percentage = percentage / 100
    count = int(len(list_a) * percentage)
    if not count: return []  # edge case, no elements removed
    list_a[-count:], list_b = [], list_a[-count:]
    return list_b


calcArgs0 = list(itertools.product(inputFileList, outPrefix, calcStat[0], modelFile, model, windowSize, windowStep, overwrite))
calcArgs1 = list(itertools.product(inputFileList, outPrefix, calcStat[1], modelFile, model, overwrite))
calcArgs2 = list(itertools.product(inputFileList, outPrefix, calcStat[2], windowSize, overwrite))
calcArgs3 = list(itertools.product(inputFileList, outPrefix, calcStat[3], windowSize, windowStep, overwrite))
calcArgs4 = list(itertools.product(inputFileList, outPrefix, calcStat[4], overwrite))
calcArgs5 = list(itertools.product(inputFileList, outPrefix, calcStat[5], overwrite))
calcArgs6 = list(itertools.product(inputFileList, outPrefix, calcStat[6], overwrite))
calcArgs7 = list(itertools.product(inputFileList, outPrefix, calcStat[7], modelFile, model, overwrite))

all_combinations = list_of_lists(calcArgs0) + list_of_lists(calcArgs1) + list_of_lists(calcArgs2) + list_of_lists(calcArgs3) + list_of_lists(calcArgs4) + list_of_lists(calcArgs5) + list_of_lists(calcArgs6) + list_of_lists(calcArgs7)


@pytest.mark.parametrize('calc_args', all_combinations, ids=id_func(all_combinations))
@pytest.mark.timeout(TIMEOUT)      # 15 minutes test case timeout
def test_calc(calc_args):
    p = Popen(['python ' + os.getcwd() + '/vcf_calc.py ' + ' '.join(flatten(calc_args))], shell=True, stdin=PIPE,
              stdout=PIPE, stderr=STDOUT, close_fds=True)

    q = psutil.Process(p.pid)
    mem = q.memory_info()[0] / 2. ** 10
    cpu = q.cpu_percent(interval=0.1)
    print('Memory Usage:', mem, 'kB')
    print('CPU Usage:', cpu, '%')

    output = p.stdout.read().decode().split('\n')
    p.communicate()
    assert p.returncode == 0, output


# pytest --html=reports/vcf_calc_test_suite.html --capture=fd --self-contained-html --capture=sys --tb=no
