# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:57:52 2017

@author: lfiorito
"""
import sys
import logging
from sandy.files import File

class Sampling(File):
    """
    Class dedicated to the sampling procedure.
    """
    
    def __repr__(self):
        return "<SAMPLING>"
    
    def __init__(self, file):
        logging.info("\n"+80*"=")
        logging.info("{:^80}".format("SAMPLING"))
        logging.info(80*"=")
        from sandy.endf import ENDF
        self.filename = file
        self.sampled_sections = []
        try:
            self.F = ENDF.read_endf(self.filename)
        except SystemExit:
            logging.error(self.prefix + "File '{}' is not processed as an error was raised while parsing".format(self.filename))
            return
        self.run()

    @property
    def prefix(self):
        return "SAMPLING : "

    @property
    def nsmp(self):
        r"""
        Number of samples.
        """
        from sandy.sandy_input import options
        if 'samples' not in options:
            logging.error(self.prefix + "Missing input keyword 'samples'")
            sys.exit()
        return options['samples']
    
    @property
    def mf(self):
        r"""
        ``MF`` sections to process.
        """
        from sandy.sandy_input import options
        return options['mf'] if 'mf' in options else self.F.keys()
    
    @property
    def sandy_file(self):
        r"""
        ``.sandy`` output file.
        """
        from sandy.sandy_input import outdir
        from os.path import join
        name = self.F.name + '.sandy'
        return join(outdir, name)

    def run(self):
        r"""
        Run sampling calculation.
        """
        from sandy.sandy_input import write
        if 31 in self.mf:
            self.sample_mf31()
        if 32 in self.mf:
            self.sample_mf32()
        if 33 in self.mf:
            self.sample_mf33()
        if 34 in self.mf:
            self.sample_mf34()
        if 35 in self.mf:
            self.sample_mf35()
        if not self.sampled_sections:
            logging.info(self.prefix + "Sampling was not performed as no requested covariance was found")
            return
        if write:
            logging.debug(80*"+")
            logging.debug("{:^80}".format(" Writing samples "))
            logging.debug(80*"+")
            self.write_sandy()
            self.write_samples()
    
    def write_samples(self):
        r"""
        Write sampled data into perturbed ``ENDF-6`` files.
        """
        from sandy.sandy_input import outdir
        from os.path import exists, join
        from os import makedirs
        from sandy.functions import printProgress
        for ismp in range(self.nsmp):
            printProgress(ismp+1, self.nsmp, prefix = 'WRITING :', suffix = 'Complete', barLength = 50)
            suffix = 'smp-{}'.format(ismp+1)
            namesmp = self.F.name + '.' + suffix
            directory = join(outdir, suffix)
            if not exists(directory): # If the directory does not exist, create it
                makedirs(directory)
            file = join(directory, namesmp)
            self.F.write_to_file(file, ismp=ismp)
    
    def write_sandy(self):
        r"""
        Write best estimate ``ENDF-6`` file after having read and processed 
        the data.
        This method is called in debug mode and can be used to check if 
        the output format is properly written.
        """
        if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
            logging.info(self.prefix + "Writing best estimate file '{}'".format(self.sandy_file))
            self.F.write_to_file(self.sandy_file)

    def sample_mf31(self):
        r"""
        Draw samples for fission neutron multiplicities.
        """
        from sandy.sandy_input import write, plot
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Sampling from MF31 "))
        logging.debug(80*"+")
        if 31 in self.F.keys():
            self.F[31].sampling()
            if write:
                self.F[31].write_to_csv(self.F.name)
            self.F[1].add_samples(self.F[31])
            if plot:
                self.F[1].plot_samples(self.F.name)
            self.sampled_sections.append(31)
        else:
            logging.warn(self.prefix + "Section 'MF31' was not found")

    def sample_mf32(self):
        r"""
        Draw samples for resonance parameters and scattering radii.
        """
        from sandy.sandy_input import write
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Sampling from MF32 "))
        logging.debug(80*"+")
        if 32 in self.F.keys():
            self.F[2][151].sample_resonances(self.F[32][151])
            self.F[2][151].sample_sr(self.F[32][151])
            if write:
                self.F[32][151].write_to_csv(self.F.name)
            self.sampled_sections.append(32)
        else:
            logging.info(self.prefix + "Section 'MF32' was not found")

    def sample_mf33(self):
        r"""
        Draw samples for cross sections.
        """
        from sandy.sandy_input import write, plot
        from sandy.endf import ENDF
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Sampling from MF33 "))
        logging.debug(80*"+")
        if 33 in self.F.keys():
            self.F[33].sampling()
            if write:
                self.F[33].write_to_csv(self.F.name)
            FP = ENDF.read_pendf(self.F)
            if 152 in FP[2]:
                logging.debug(self.prefix + "Copy section MF2/MT152 from PENDF")
                self.F[2].update({ 152 : FP[2][152]})
            logging.debug(self.prefix + "Copy section MF3 from PENDF")
            self.F[3] = FP[3]
            self.F[3].add_samples(self.F[33])
            if plot:
                self.F[3].plot_samples(self.F.name)
            self.F.INFO.lrp = 2
            self.sampled_sections.append(33)
        else:
            logging.warn(self.prefix + "Section 'MF33' was not found")

    def sample_mf34(self):
        r"""
        Draw samples for angular distributions.
        """
        from sandy.sandy_input import write, plot
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Sampling from MF34 "))
        logging.debug(80*"+")
        if 34 in self.F.keys():
            self.F[34].sampling()
            if write:
                self.F[34].write_to_csv(self.F.name)
            self.F[4].add_samples(self.F[34])
            if plot:
                self.F[4].plot_samples(self.F.name)
            self.sampled_sections.append(34)
        else:
            logging.warn(self.prefix + "Section 'MF34' was not found")

    def sample_mf35(self):
        r"""
        Draw samples for energy distributions.
        """
        from sandy.sandy_input import write, plot
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Sampling from MF35 "))
        logging.debug(80*"+")
        if 35 in self.F.keys():
            self.F[35].sampling()
            if write:
                self.F[35].write_to_csv(self.F.name)
            self.F[5].add_samples(self.F[35])
            if plot:
                self.F[5].plot_samples(self.F.name)
            self.sampled_sections.append(35)
        else:
            logging.warn(self.prefix + "Section 'MF35' was not found")



    

def run(file):
    from sandy.sandy_input import outdir, plot, write, options
    from sandy.functions import printProgress
    from sandy.endf import ENDF
    from os.path import join, split, exists
    from os import makedirs

    # CHECK KEYWORDS
    if 'samples' not in options:
        logging.error("ERROR : missing input keyword 'samples'")
        sys.exit()
    nsmp = options['samples']

    # PROCESS FILE
    path, name = split(file)
    try:
        F = ENDF.read_endf(file, **options)
    except SystemExit:
        logging.error("ERROR : File '{}' is not processed as an error was raised while parsing".format(file))
        return
    
    mf = options['mf'] if 'mf' in options else F.keys()

    found_samples= []
    if 31 in mf:
        if 31 in F.keys():
            found_samples.append(31)
            F[31].sampling(nsmp, **options)
            if 'stdmax' in options:
                F[31].delete_samples_above_threshold(options['stdmax'])
            if write:
                F[31].write_to_csv(join(outdir, F.name), **options)
            F[31].apportion_to_daughters()
            F[1].add_samples(F[31], **options)
            if plot:
                F[1].plot_samples(F.name)
        else:
            logging.warn("MF31: File 'MF31' was not found")
    if 32 in mf:
        if 32 in F.keys():
            found_samples.append(32)
            F[2][151].sample_resonances(F[32][151])
            F[2][151].sample_sr(F[32][151])
            if write:
                F[32][151].write_to_csv(F.name)
        else:
            logging.warn("MF32: File 'MF32' was not found")
    if 33 in mf:
        if 33 in F.keys():
            found_samples.append(33)
            F[33].sampling(nsmp, **options)
            if 'stdmax' in options:
                F[33].delete_samples_above_threshold(options['stdmax'])
            if write:
                F[33].write_to_csv(join(outdir, F.name), **options)
            F[33].apportion_to_daughters()
            FP = ENDF.read_pendf(F, cwd=outdir, **options)
            FP.process(**options)
            F[3] = FP[3]
            del FP
            F[3].add_samples(F[33], **options)
            if plot:
                F[3].plot_samples(F.name)
            F.INFO.lrp = 2
        else:
            logging.warn("MF33: File 'MF33' was not found")
    if 34 in mf:
        if 34 in F.keys():
            found_samples.append(34)
            F[34].sampling(nsmp, **options)
            if 'stdmax' in options:
                F[34].delete_samples_above_threshold(options['stdmax'])
            if write:
                F[34].write_to_csv(join(outdir, F.name), **options)
            F[4].add_samples(F[34], **options)
            if plot:
                F[4].plot_samples(outdir, F.name, **options)
        else:
            logging.warn("MF34: File 'MF34' was not found")
    if 35 in mf:
        if 35 in F.keys():
            found_samples.append(35)
            F[35].sampling(nsmp, **options)
            if write:
                F[35].write_to_csv(join(outdir, F.name), **options)
            F[5].add_samples(F[35], **options)
            if 'stdmax' in options:
                F[5].delete_samples_above_threshold(options['stdmax'])
            if plot:
                F[5].plot_samples(F.name)
        else:
            logging.warn("MF35: File 'MF35' was not found")

    sname = F.name + '.sandy'
    logging.info("SAMPLING : Writing best estimate file '{}'".format(join(outdir,sname)))
    F.write_to_file(sname, cwd=outdir)

    if not found_samples:
        logging.info("SAMPLING : sampling was not performed as no requested covariance was found")
        return


    # WRITE SAMPLES
    if write:
        for ismp in range(nsmp):
            printProgress(ismp+1, nsmp, prefix = 'Writing:', suffix = 'Complete', barLength = 50)
            suffix = 'smp-' + str(ismp+1)
            namesmp = F.name + '.' + suffix
            directory = join(outdir, suffix)
            if not exists(directory): # If the directory does not exist, create it
                makedirs(directory)
            F.write_to_file(name=namesmp, cwd=directory, ismp=ismp)



def preprocessing():
    r"""
    Preprocess `SANDY`'s inputs before running in sampling mode.
    """
    from sandy.sandy_input import process_input, get_endf_files
    inp = process_input()
    files = get_endf_files(inp)
    for file in files:
        logging.info("SAMPLING : Run sampling module for file '{}'".format(file))
        run(file, **inp)
    logging.info("CALCULATION COMPLETED")



if __name__ == "__main__":
    preprocessing()