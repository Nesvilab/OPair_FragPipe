using System;
using System.Collections.Generic;
using System.IO;

namespace EngineLayer
{
    public class FilterRule
    {
        public FilterRule(List<double> massList, Dictionary<byte, int> rules) 
        { 
            Masses = massList;
            MonosaccharidesRequired = rules;
        }

        public List<double> Masses { get; }

        public Dictionary<byte, int> MonosaccharidesRequired { get; }

        /**
        * Load a list of rules for oxonium ion-based filtering from a tsv file. 
        * File format: 
        * mass1,mass2,...\t ID1 \t minCount1 \t ID2 \t minCount2 (etc)
        */
        public static List<FilterRule> LoadOxoniumFilters(string filePath)
        {
            List<FilterRule> rules = new();
            using (var lines = new StreamReader(filePath))
            {
                var firstLine = true;
                while (lines.Peek() != -1)
                {
                    var line = lines.ReadLine();

                    if (firstLine)
                    {
                        //Skip the first line
                        firstLine = false;
                        continue;
                    }

                    var splits = line.Replace('"', ' ').Split('\t');
                    
                    // for backwards compatibility
                    var startCol = 1;
                    if (splits[0].Contains('('))
                    {
                        startCol = 0;   // new version has one less column
                    }

                    // parse rules from required residue composition
                    var kind = GlycanDatabase.String2Kind(splits[startCol]);
                    Dictionary<byte, int> ruleDict = new();
                    for (var i = 0; i < kind.Length; i++)
                    {
                        if (kind[i] > 0)
                        {
                            ruleDict[(byte) i] = kind[i];
                        }
                    }
                    
                    // Parse mass
                    List<double> masses = new();
                    var glycan = new Glycan(kind);
                    var startMass = glycan.Mass / 1E5 + Chemistry.Constants.ProtonMass;
                    startMass += Double.Parse(splits[startCol + 1]);
                    masses.Add(startMass);

                    var rule = new FilterRule(masses, ruleDict);
                    rules.Add(rule);
                }
            }
            return rules;
        }
    }
}
