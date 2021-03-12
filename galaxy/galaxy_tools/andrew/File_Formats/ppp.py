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

# Currently these supported binary data types must be manually set on upload

class VcfGz(Binary):
    """Class describing a VCF.GZ file"""
    #edam_format = ""
    #edam_data = ""
    file_ext = "vcf.gz"

    compressed = True

    MetadataElement( name="vcfgz_index", desc="VCF.GZ Index File", param=metadata.FileParameter, file_ext="tbi", readonly=True, no_value=None, visible=False, optional=True )

    def sniff( self, filename ):
        # VCF.GZ is compressed in the BGZF format, and must not be uncompressed in Galaxy.
        try:
            f = open(vcfname,'rb')
            l = f.readline()
            f.close()
            BGZF_HEADER=b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00'
            if l[:len(BGZF_HEADER)] == BGZF_HEADER:
                return True
            return False
        except:
            return False

    def set_meta( self, dataset, overwrite=True, **kwd ):
        """ Creates the index for the VCF.GZ file. """
        # These metadata values are not accessible by users, always overwrite
        index_file = dataset.metadata.vcfgz_index
        if not index_file:
            index_file = dataset.metadata.spec['vcfgz_index'].param.new_file( dataset=dataset )
        # Create the vcf.gz index
        # $ bcftools index
        # Usage: bcftools index -t <in.vcf.gz>

        dataset_symlink = os.path.join( os.path.dirname( index_file.file_name ),
                                        '__dataset_%d_%s' % ( dataset.id, os.path.basename( index_file.file_name ) ) )
        os.symlink( dataset.file_name, dataset_symlink )

        stderr_name = tempfile.NamedTemporaryFile( prefix="vcfgz_index_stderr" ).name
        command = [ 'bcftools', 'index', '-t', dataset_symlink ]
        try:
            subprocess.check_call( args=command, stderr=open( stderr_name, 'wb' ) )
            shutil.move( dataset_symlink + '.tbi', index_file.file_name )
        except Exception as e:
            stderr = open( stderr_name ).read().strip()
            raise Exception('Error setting VCF.GZ metadata: %s' % (stderr or str(e)))
        finally:
            # Remove temp file and symlink
            os.remove( stderr_name )
            os.remove( dataset_symlink )
        dataset.metadata.vcfgz_index = index_file


Binary.register_sniffable_binary_format("vcf.gz", "vcf.gz", VcfGz)


class Model(Json):
    file_ext = ".model"
    #edam_format = "format_3746"

    MetadataElement(name="models", default=[], desc="Models", param=MetadataParameter, readonly=True, visible=False, optional=True, no_value=[])

    # Update as format becomes more defined
    def set_peek(self, dataset, is_multi_byte=False):
        if not dataset.dataset.purged:
            dataset.blurb = "Model File"

    # Update as format becomes more defined
    def sniff(self, filename):
        if self._looks_like_json(filename):
            return True
        else:
            return False

    def set_meta(self, dataset, **kwd):
        if dataset.has_data():
            with open(dataset.file_name) as model_file:
                try:
                    model_dict = json.load(model_file)
                except Exception:
                    return

                models = [model['name'] for model in model_dict if 'name' in model]
        
                dataset.metadata.models = models
