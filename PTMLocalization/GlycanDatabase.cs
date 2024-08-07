﻿using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;

namespace EngineLayer
{

    public static class GlycanDatabase
    {
        // TODO: read from monosaccharides file in case it is changed!
        //kinds = Monosaccharide.GetCharMassDic(GlobalVariables.Monosaccharides).Keys.ToString();
        private static string kinds;
        private static Regex SimpleKindPattern = new Regex(@"[NHAGFPSYCXUM][0-9]+");

        /**
         * Handler for file or glycan list passed to this parameter. Detect if file or glycan list and parse accordingly.
         */
        public static IEnumerable<Glycan> LoadGlycanDatabase(string fileOrGlycans, bool ToGenerateIons, bool IsLoadOGlycan)
        {
            if (File.Exists(fileOrGlycans))
            {
                return LoadGlycan(fileOrGlycans, ToGenerateIons, IsLoadOGlycan);
            }
            else
            {
                // the provided string is not a file, so load it as a glycan list
                return LoadGlycanList(fileOrGlycans, ToGenerateIons, IsLoadOGlycan);
            }
        }

        //Load Glycan. Generally, glycan-ions should be generated for N-Glycopepitdes which produce Y-ions; MS method couldn't produce o-glycan-ions.
        public static IEnumerable<Glycan> LoadGlycan(string filePath, bool ToGenerateIons, bool IsLoadOGlycan)
        {
            bool isKind = true;
            bool isSimpleKind = false;
            using (StreamReader lines = new StreamReader(filePath))
            {
                while(lines.Peek() != -1)
                {
                    string line = lines.ReadLine();
                    // skip headers
                    if (line.StartsWith("#") || line.StartsWith("//") || line.StartsWith("%") || line.StartsWith("Glycan") || line.Contains(','))
                    {
                        line = lines.ReadLine();
                    }
                    if (!line.Contains("HexNAc"))
                    {
                        // Distinguish simple composition (e.g, N1H1A1) from pGlyco style stucts (e.g., (N(H(A)) )
                        if (!line.Contains('('))
                        {
                            isSimpleKind = true;
                        } else
                        {
                            isKind = false;
                        }
                    }
                    break;
                }
            }

            if (isKind)
            {
                if (isSimpleKind)
                {
                    return LoadSimpleKindGlycan(filePath, ToGenerateIons, IsLoadOGlycan);
                }
                else
                {
                    return LoadKindGlycan(filePath, ToGenerateIons, IsLoadOGlycan);
                }
            }
            else
            {
                return LoadStructureGlycan(filePath, ToGenerateIons, IsLoadOGlycan);
            }
        }

        public static IEnumerable<Glycan> LoadGlycanList(string glycanList, bool ToGenerateIons, bool IsOGlycanSearch)
        {
            string[] splits = Regex.Split(glycanList, "[,;\\s]");
            int id = 1;
            foreach (string glycanStr in splits)
            {
                if (glycanStr.IsNotNullOrEmptyOrWhiteSpace())
                {
                    var kind = String2Kind(glycanStr);
                    var glycan = new Glycan(kind);
                    glycan.GlyId = id++;
                    if (ToGenerateIons)
                    {
                        if (IsOGlycanSearch)
                        {
                            glycan.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        else
                        {
                            glycan.Ions = NGlycanCompositionFragments(kind);
                        }
                    }

                    yield return glycan;
                }
            }
        }

        //Load KindGlycan. Compatible with Byonic.
        public static IEnumerable<Glycan> LoadKindGlycan(string filePath, bool ToGenerateIons, bool IsOGlycanSearch)
        {
            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine().Split('\t').First();

                    if (!(line.Contains("HexNAc") || line.Contains("Hex")))
                    {
                        continue;
                    }

                    var kind = String2Kind(line);

                    var glycan = new Glycan(kind);
                    glycan.GlyId = id++;
                    if (ToGenerateIons)
                    {
                        if (IsOGlycanSearch)
                        {
                            glycan.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        else
                        {
                            glycan.Ions = NGlycanCompositionFragments(kind);
                        }
                    }
                    yield return glycan;
                }
            }
        }

        //Load simple composition string. Example format: N1H1A1
        public static IEnumerable<Glycan> LoadSimpleKindGlycan(string filePath, bool ToGenerateIons, bool IsOGlycanSearch)
        {
            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine().Split('\t').First();

                    if (line.Contains('#') || line.Contains("//"))
                    {
                        continue;
                    }

                    var kind = SimpleString2Kind(line);

                    var glycan = new Glycan(kind);
                    glycan.GlyId = id++;
                    if (ToGenerateIons)
                    {
                        if (IsOGlycanSearch)
                        {
                            glycan.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        else
                        {
                            glycan.Ions = NGlycanCompositionFragments(kind);
                        }
                    }
                    yield return glycan;
                }
            }
        }

        public static byte[] String2Kind(string line)
        {
            byte[] kind = new byte[Glycan.SugarLength];
            var x = line.Split(new char[] { '(', ')' });
            int i = 0;
            while (i < x.Length - 1)
            {
                kind[Glycan.NameIdDic[x[i]]] = byte.Parse(x[i + 1]);
                i = i + 2;
            }

            return kind;
        }

        public static byte[] SimpleString2Kind(string line)
        {
            byte[] kind = new byte[Glycan.SugarLength];
            Match matcher = SimpleKindPattern.Match(line);
            while (matcher.Success)
            {
                string glycanToken = matcher.Value;
                char glycanType = Char.Parse(glycanToken.Substring(0, 1));
                int glycanCount = int.Parse(glycanToken.Substring(1));
                if (Monosaccharide.GetSymbolIdDic(GlobalVariables.Monosaccharides).ContainsKey(glycanType))
                {
                     kind[Monosaccharide.GetSymbolIdDic(GlobalVariables.Monosaccharides)[glycanType]] = (byte)glycanCount;
                }
                else
                {
                    Console.WriteLine($"Invalid token {glycanToken} in line {line}. This line will be skipped");
                    return new byte[Glycan.SugarLength];
                }
                matcher = matcher.NextMatch();
            }
            return kind;
        }

        //Load structured Glycan database.
        public static IEnumerable<Glycan> LoadStructureGlycan(string filePath, bool ToGenerateIons, bool IsOGlycan)
        {
            using (StreamReader glycans = new StreamReader(filePath))
            {
                int id = 1;
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();
                    if (line.Contains(','))
                    {
                        continue;
                    }
                    yield return Glycan.Struct2Glycan(line, id++, ToGenerateIons, IsOGlycan);
                }
            }
        }

        //This function build fragments based on the general core of NGlyco fragments. 
        //From https://github.com/mobiusklein/glycopeptidepy/structure/fragmentation_strategy/glycan.py#L408
        //The fragment generation is not as good as structure based method. So it is better to use a structure based N-Glycan database.
        public static List<GlycanIon> NGlycanCompositionFragments(byte[] kind)
        {
            int glycan_mass = Glycan.GetMass(kind);

            int core_count = 1;
            int iteration_count = 0;
            bool extended = true;
            bool extended_fucosylation = false;

            int fuc_count = kind[4];
            int xyl_count = kind[9];
            int hexnac_inaggregate = kind[0];
            int hexose_inaggregate = kind[1];

            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int base_hexnac = Math.Min(hexnac_inaggregate + 1, 3);
            for (int hexnac_count = 0; hexnac_count < base_hexnac; hexnac_count++)
            {
                if (hexnac_count == 0)
                {
                    byte[] aIonKind = new byte[Glycan.SugarLength];
                    aIonKind[1] = (byte)hexnac_count;
                    GlycanIon glycanIon = new GlycanIon(null, 8303819, aIonKind, glycan_mass - 8303819);
                    glycanIons.Add(glycanIon);
                }
                else if (hexnac_count == 1)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);

                    glycanIons.Add(glycanIon);

                    if (iteration_count < fuc_count)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);

                        glycanIons.Add(fuc_glycanIon);
                    }
                }
                else if (hexnac_count == 2)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);
                    glycanIons.Add(glycanIon);

                    if (!extended_fucosylation)
                    {
                        if (iteration_count < fuc_count)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }
                    }
                    else if (fuc_count > 0)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                        glycanIons.Add(fuc_glycanIon);

                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                        {
                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                            glycanIons.Add(add_fuc_glycanIon);
                        }

                        if (iteration_count < xyl_count)
                        {
                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                            glycanIons.Add(xyl_fuc_glycanIon);
                        }
                    }

                    if (iteration_count < xyl_count)
                    {
                        GlycanIon xyl_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 0, 1, glycan_mass);
                        glycanIons.Add(xyl_glycanIon);
                    }


                    int min_hexose_inaggregate = Math.Min(hexose_inaggregate + 1, 4);
                    for (int hexose_count = 1; hexose_count <= min_hexose_inaggregate; hexose_count++)
                    {
                        GlycanIon hexose_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)hexnac_count, 0, 0, glycan_mass);
                        glycanIons.Add(hexose_glycanIon);

                        if (!extended_fucosylation)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }
                        else if (fuc_count > 0)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                            {
                                GlycanIon add_fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                                glycanIons.Add(add_fuc_glycanIon);
                            }

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }

                        if (iteration_count < xyl_count)
                        {
                            GlycanIon xyl_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 0, 1, glycan_mass);
                            glycanIons.Add(xyl_glycanIon);
                        }

                        if (hexose_count == 3 && hexnac_count >= 2 * core_count && extended)
                        {
                            for (int extra_hexnac_count = 0; extra_hexnac_count < hexnac_inaggregate - hexnac_count + 1; extra_hexnac_count++)
                            {
                                if (extra_hexnac_count + hexnac_count > hexnac_inaggregate)
                                {
                                    continue;
                                }

                                if (extra_hexnac_count > 0)
                                {
                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);

                                    if (!extended_fucosylation)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }
                                    else if (fuc_count > 0)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                                        {
                                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                                            glycanIons.Add(add_fuc_glycanIon);
                                        }

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }

                                    if (iteration_count < xyl_count)
                                    {
                                        GlycanIon xyl_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 0, 1, glycan_mass);
                                        glycanIons.Add(xyl_glycanIon);
                                    }

                                }

                                for (int extra_hexose_count = 1; extra_hexose_count < hexose_inaggregate - hexose_count + 1; extra_hexose_count++)
                                {
                                    if (extra_hexose_count + hexose_count > hexose_inaggregate)
                                    {
                                        continue;
                                    }

                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)(hexose_count + extra_hexose_count), (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);

                                    if (!extended_fucosylation)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }
                                    else if (fuc_count > 0)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                                        {
                                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                                            glycanIons.Add(add_fuc_glycanIon);
                                        }

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }

                                    if (iteration_count < xyl_count)
                                    {
                                        GlycanIon xyl_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 0, 1, glycan_mass);
                                        glycanIons.Add(xyl_glycanIon);
                                    }

                                }
                            }
                        }
                    }

                }


            }

            return glycanIons;
        }

        private static GlycanIon GenerateGlycanIon(byte hexose_count, byte hexnac_count, byte fuc_count, byte xyl_count, int glycan_mass)
        {
            byte[] aIonKind = new byte[Glycan.SugarLength];
            aIonKind[0] = hexose_count; aIonKind[1] = hexnac_count; aIonKind[4] = fuc_count; aIonKind[9] = xyl_count;

            int ionMass = Glycan.GetMass(aIonKind);

            GlycanIon glycanIon = new GlycanIon(null, ionMass, aIonKind, glycan_mass - ionMass);

            return glycanIon;
        }

        private static GlycanIon ExtendGlycanIon(GlycanIon glycanIon, byte hexose_count, byte hexnac_count, byte fuc_count, byte xyl_count, int glycan_mass)
        {
            byte[] ionKind = glycanIon.IonKind;
            ionKind[0] += hexose_count;
            ionKind[1] += hexnac_count;
            ionKind[4] += fuc_count;
            ionKind[9] += xyl_count;

            int ionMass = Glycan.GetMass(ionKind);

            GlycanIon extend_glycanIon = new GlycanIon(null, ionMass, ionKind, glycan_mass - ionMass);

            return extend_glycanIon;
        }

        //This function build fragments based on the general core of OGlyco fragments. 
        //From https://github.com/mobiusklein/glycopeptidepy/structure/fragmentation_strategy/glycan.py
        //The fragment generation is not as good as structure based method. So it is better to use a structure based O-Glycan database.
        public static List<GlycanIon> OGlycanCompositionFragments(byte[] kind)
        {
            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int glycan_mass = Glycan.GetMass(kind);

            int iteration_count = 0;
            bool extended = true;

            int fuc_count = kind[4];
            int hexnac_inaggregate = kind[0];
            int hexose_inaggregate = kind[1];

            for (int hexnac_count = 0; hexnac_count < 3; hexnac_count++)
            {
                if (hexnac_inaggregate < hexnac_count)
                {
                    continue;
                }


                if (hexnac_count >= 1)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);

                    glycanIons.Add(glycanIon);

                    if (iteration_count < fuc_count)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);

                        glycanIons.Add(fuc_glycanIon);
                    }

                    for (int hexose_count = 0; hexose_count < 2; hexose_count++)
                    {
                        if (hexose_inaggregate < hexose_count)
                        {
                            continue;
                        }

                        if (hexose_count > 0)
                        {
                            GlycanIon hexose_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)hexnac_count, 0, 0, glycan_mass);
                            glycanIons.Add(hexose_glycanIon);

                            if (iteration_count < fuc_count)
                            {
                                GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);

                                glycanIons.Add(fuc_glycanIon);
                            }
                        }

                        // After the core motif has been exhausted, speculatively add on the remaining core monosaccharides sequentially until exhausted.

                        if (extended && hexnac_inaggregate - hexnac_count >= 0)
                        {
                            for (int extra_hexnac_count = 0; extra_hexnac_count  < hexnac_inaggregate - hexnac_count + 1; extra_hexnac_count ++)
                            {
                                if (extra_hexnac_count > 0)
                                {
                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);


                                    if (iteration_count < fuc_count)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);

                                        glycanIons.Add(fuc_glycanIon);
                                    }

                                }

                                if (hexose_inaggregate > hexose_count && hexose_count > 0)
                                {
                                    for (int extra_hexose_count = 0; extra_hexose_count < hexose_inaggregate - hexose_count; extra_hexose_count++)
                                    {
                                        if (extra_hexose_count > 0 && extra_hexose_count + hexose_count >0)
                                        {

                                            GlycanIon new_glycanIon = GenerateGlycanIon((byte)(hexose_count + extra_hexose_count), (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                            glycanIons.Add(new_glycanIon);


                                            if (iteration_count < fuc_count)
                                            {
                                                GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);

                                                glycanIons.Add(fuc_glycanIon);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }

            }


            return glycanIons;
        }

        //The OGlycanCompositionFragments just generate some core GlycanIons. We need a combination solution.
        public static List<GlycanIon> OGlycanCompositionCombinationChildIons(byte[] kind)
        {
            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int glycan_mass = Glycan.GetMass(kind);

            List<byte[]> _kinds = new List<byte[]>();
            HashSet<string> _keys = new HashSet<string>();

            _kinds.Add((byte[])kind.Clone());
            _GetCombinations(kind, _kinds, _keys);

            foreach (var k in _kinds)
            {
                //Rules to build OGlycan child ions.
                //At least one HexNAc
                if (k[1] == 0)
                {
                    continue;
                }

                //#Fucose <= #HexNAc. One Fucose modify one 
                if (k[4]!= 0 && k[4] > k[1] )
                {
                    continue;
                }

                //#NeuAc * 2 >= #Acetylation. One NeuAc can be modified with two Acetylation
                if (k[9]!= 0 && k[2]*2 < k[9])
                {
                    continue;
                }

                var ionMass = Glycan.GetMass(k);
                GlycanIon glycanIon = new GlycanIon(null, ionMass, k, glycan_mass - ionMass);
                glycanIons.Add(glycanIon);
            }

            return glycanIons.OrderBy(p=>p.IonMass).ToList();
        }

        private static void _GetCombinations(byte[] kind, List<byte[]> _kinds, HashSet<string> _keys)
        {
            if (kind.Sum(p=>p) == 0)
            {
                return;
            }
            else
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    if (kind[i] >= 1)
                    {
                        byte[] akind = (byte[])kind.Clone();
                        akind[i]--;
                        if (akind.Sum(p => p) != 0)
                        {
                            if (!_keys.Contains(Glycan.GetKindString(akind)))
                            {
                                _keys.Add(Glycan.GetKindString(akind));
                                _kinds.Add((byte[])akind.Clone());
                                _GetCombinations(akind, _kinds, _keys);
                            }
                           
                        }

                    }
                }
            }
        }
    }
}
