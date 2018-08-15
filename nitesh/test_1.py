import os
import pytest
import psutil
import itertools
from subprocess import Popen, PIPE, STDOUT

TIMEOUT = 60 * 15  # 15 minute timeout holder

inputFileList = [['--vcf', 'chr11.unfiltered.vcf']]
outFormatList = [['--out-format', 'vcf'], ['--out-format', 'vcf.gz'], ['--out-format', 'bcf'], ['--out-format', 'removed_sites'], ['--out-format', 'kept_sites'], ['--out-format', 'removed_bed'], ['--out-format', 'kept_bed']]
outPrefixList = [['--out-prefix', 'chr11_vcf']]
Remaining_options1 = [None, ['--filter-include-indv', 'Pan_paniscus-9731_LB502'], ['--filter-exclude-indv', 'Pan_paniscus-A914_Hortense'], ['--filter-include-indv-file', 'Paniscus.txt'], ['--filter-exclude-indv-file', 'Schweinfurthii.txt'], ['--filter-min-alleles', '2'], ['--filter-max-alleles', '4'], ['--filter-max-missing', '1.0'], ['--filter-include-chr', 'chr1'], ['--filter-exclude-chr', 'chr2']]
Remaining_options2 = [None, ['--filter-from-bp', '1000000'], ['--filter-to-bp', '2000000'], ['--filter-include-positions', 'chr1_sites.sites'], ['--filter-exclude-positions', 'chr1_sites.sites'], ['--filter-include-bed', 'chr1_sites.bed'], ['--filter-exclude-bed', 'chr1_sites.bed'], ['--filter-include-passed'], ['--filter-include-flag', 'q10']]
Remaining_options3 = [None, ['--filter-exclude-flag', 's50'], ['--filter-include-info', 'DB'], ['--filter-exclude-info', 'H2'], ['--filter-maf-min', '0.2'], ['--filter-maf-max', '0.3'], ['--filter-mac-min', '8'], ['--filter-mac-max', '12'], ['--filter-distance', '500']]
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
        ids.append('python vcf_filter.py ' + ' '.join(i))
    return ids


filterArgs0 = list(itertools.product(inputFileList, outFormatList, outPrefixList, Remaining_options1, Remaining_options2, Remaining_options3, overwrite))
filterArgs1 = list(itertools.product(inputFileList, outFormatList, outPrefixList, list_of_lists(itertools.combinations(filter(None.__ne__, Remaining_options1), 2)), overwrite))
filterArgs2 = list(itertools.product(inputFileList, outFormatList, outPrefixList, list_of_lists(itertools.combinations(filter(None.__ne__, Remaining_options2), 2)), overwrite))
filterArgs3 = list(itertools.product(inputFileList, outFormatList, outPrefixList, list_of_lists(itertools.combinations(filter(None.__ne__, Remaining_options3), 2)), overwrite))

all_combinations = list_of_lists(filterArgs0) + list_of_lists(filterArgs1) + list_of_lists(filterArgs2) + list_of_lists(filterArgs3)


@pytest.mark.parametrize('filter_args', all_combinations, ids=id_func(all_combinations))
@pytest.mark.timeout(TIMEOUT)      # 15 minutes test case timeout
def test_filter(filter_args):
    p = Popen(['python ' + os.getcwd() + '/vcf_filter.py ' + ' '.join(flatten(filter_args))], shell=True, stdin=PIPE,
              stdout=PIPE, stderr=STDOUT, close_fds=True)

    q = psutil.Process(p.pid)
    mem = q.memory_info()[0] / 2. ** 10
    cpu = q.cpu_percent(interval=0.1)
    print('Memory Usage:', mem, 'kB')
    print('CPU Usage:', cpu, '%')

    output = p.stdout.read().decode().split('\n')
    p.communicate()
    assert p.returncode == 0, output


# pytest.main(['-n', 'auto', '--html=reports/vcf_filter_test_suite.html', '--capture=fd', '--self-contained-html', '--capture=sys', '--tb=no'])
# pytest -n auto --html=reports/vcf_filter_test_suite.html --capture=fd --self-contained-html --capture=sys --tb=no
