﻿using System;
using System.Linq;
using MzLibUtil;
using EngineLayer;
using System.IO;

namespace PTMLocalization
{
    public class Task
    {
        public int run_msfragger(double productPpmTol, double PrecursorPpmTol, string psmFile, string scanpairFile, 
            string rawfileDirectory, string lcmsFileList, string glycoDatabase, int maxNumGlycans, int minIsotopeError, 
            int maxIsotopeError, bool filterOxonium, double oxoMinRelativeIntensity, int numThreads, string[] allowedSites)
        {
            Tolerance ProductMassTolerance = new PpmTolerance(productPpmTol);
            Tolerance PrecursorMassTolerance = new PpmTolerance(PrecursorPpmTol);
            if (psmFile == null)
            {
                Console.WriteLine("No PSM file specified, exiting.");
                return 2;
            }


            if (!File.Exists(glycoDatabase))
            {
                // read internal databases if not passed a full path
                try
                {
                    glycoDatabase = GlobalVariables.OGlycanLocations.Where(p => p.Contains(glycoDatabase)).First();
                }
                catch (Exception ex)
                {
                    // ignore - glycans passed as list, not file
                }

            }

            int[] isotopes = new int[maxIsotopeError - minIsotopeError + 1];
            for (int i = 0; i < isotopes.Length; i++)
            {
                isotopes[i] = minIsotopeError + i;
            }

            var localizer = new MSFragger_RunLocalization(psmFile, scanpairFile, rawfileDirectory, lcmsFileList, glycoDatabase, 
                maxNumGlycans, PrecursorMassTolerance, ProductMassTolerance, isotopes, filterOxonium, oxoMinRelativeIntensity, numThreads, allowedSites);
            int returnCode = localizer.Localize();
            return returnCode;
        }

    }
}
