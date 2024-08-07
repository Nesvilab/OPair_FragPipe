﻿using CommandLine;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;


namespace CMD
{
    public class CmdSettings
    {

        [Option('v', Default = VerbosityType.normal, HelpText = "[Optional] Determines how much text is written. Options are no output ('none'), minimal output and errors  ('minimal'), or normal ('normal')")]
        public VerbosityType Verbosity { get; set; }
        public enum VerbosityType { none, minimal, normal };

        [Option('b', Default = 10, HelpText = "[Optional] productPpmTol ")]
        public double productPpmTol { get; set; }

        [Option('c', Default = 30, HelpText = "[Optional] precursorPpmTol ")]
        public double precursorPpmTol { get; set; }

        [Option('g', Default = "OGlycan.gdb", HelpText = "[Optional] glycoDatabase ")]
        public string glycoDatabase { get; set; }

        [Option('x', Default = "", HelpText = "[Optional] glycanResiduesDatabase ")]
        public string glycanResiduesFile { get; set; }

        [Option('y', Default = "", HelpText = "[Optional] glycanModsDatabase ")]
        public string glycanModsFile { get; set; }

        [Option('n', Default = 3, HelpText = "[Optional] maxNumGlycans ")]
        public int maxNumGlycans { get; set; }

        [Option('i', Default = 0, HelpText = "[Optional] minIsotopeError ")]
        public int minIsotopeError { get; set; }

        [Option('j', Default = 2, HelpText = "[Optional] maxIsotopeError ")]
        public int maxIsotopeError { get; set; }

        [Option('d', HelpText = "[Optional] rawfileDirectory")]
        public string rawfileDirectory { get; set; }

        [Option('r', HelpText = "[Optional] LCMS file list (supercedes rawfileDirectory (-d))")]
        public string lcmsFilesList { get; set; }

        [Option('s', HelpText = "[Optional] psmFile")]
        public string psmFile { get; set; }

        [Option('p', HelpText = "[Optional] scanpairFile. Searches raw file directory if not provided ")]
        public string scanpairFile { get; set; }

        [Option('o', HelpText = "[Optional] Output folder")]
        public string outputFolder { get; set; }

        [Option('f', HelpText = "[Optional] If present, filter glycans based on oxonium ions observed. Uses default list if path to a tsv file is not provided")]
        public string oxoFilter { get; set; }

        [Option('m', Default = 0, HelpText = "[Optional] Minimum relative intensity for oxonium ion filtering. Summed intensity for all oxonium ions in a class must exceed this ratio relative to the spectrum base peak to pass filtering. Values between 0 and 1.")]
        public double oxoMinInt { get; set; }

        [Option('t', Default = 0, HelpText = "[Optional] Number of threads")]
        public int numThreads { get; set; }

        public override string ToString()
        {
            string x = "";
            x += " -t " + numThreads.ToString();
            x += " -b " + productPpmTol.ToString();

            x += " -c " + precursorPpmTol.ToString();

            x += " -g " + glycoDatabase;
            x += " -x " + glycanResiduesFile;
            x += " -y " + glycanModsFile;
            x += " -n " + maxNumGlycans.ToString();
            x += " -i " + minIsotopeError.ToString(); 
            x += " -j " + maxIsotopeError.ToString();
            x += " -f " + oxoFilter;
            x += " -m " + oxoMinInt.ToString();
            x += " -d " + rawfileDirectory;
            x += " -r " + lcmsFilesList;
            x += " -s " + psmFile;
            x += " -p " + scanpairFile;
            x += " -o " + outputFolder;

            return x;
        }
    }
}
