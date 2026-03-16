

namespace EngineLayer.GlycoSearch
{
    // Localization of multiple glycans on one peptide can be divided into the following groups based on the quality of the localization. Similar to Proteoform Level.
    public enum LocalizationLevel
    {
        Level1,
        Level1b,
        Level2,
        Level3
    }

    public enum GlycoType
    {
        SinglePep,
        OGlycoPep,
        NGlycoPep,
        MixedGlycoPep
    }

    public enum GlycoSearchType
    {
        OGlycanSearch,
        NGlycanSearch,
        N_O_GlycanSearch

    }
}
