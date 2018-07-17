from __future__ import print_function

import gzip
import logging
import os
import shutil
import subprocess
import tempfile
import json

import pysam
from bx.seq.twobit import TWOBIT_MAGIC_NUMBER, TWOBIT_MAGIC_NUMBER_SWAP, TWOBIT_MAGIC_SIZE

from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import DictParameter, ListParameter, MetadataElement, MetadataParameter

from binary import Binary
from text import Json

log = logging.getLogger(__name__)

class Model(Json):
    file_ext = ".model"
    #edam_format = "format_3746"

    MetadataElement(name="models", default=[], desc="Models", param=MetadataParameter, readonly=True, visible=False, optional=True, no_value=[])
    MetadataElement(name="npop_dict", default={}, desc="Model population count", param=MetadataParameter, readonly=True, visible=False, optional=True, no_value={})

    # Update as format becomes more defined
    def set_peek(self, dataset, is_multi_byte=False):
        if not dataset.dataset.purged:
            dataset.blurb = "Model File"

    # Update as format becomes more defined
    '''
    def sniff(self, filename):
        try:
            json.loads(filename.file_name)
            return True
        except Exception:
            return False
    '''

    def set_meta(self, dataset, **kwd):
        if dataset.has_data():
            with open(dataset.file_name) as model_file:
                try:
                    model_list = json.load(model_file)
                except Exception:
                    return

                # Assign models
                models = [model_dict['name'] for model_dict in model_list if 'name' in model_dict]

                # Assign population count per model
                npop_dict = dict([(model_dict['name'], len(model_dict['pops'])) for model_dict in model_list if 'name' in model_dict])

                dataset.metadata.models = models
                dataset.metadata.npop_dict = npop_dict
