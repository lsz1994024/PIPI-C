package proteomics;

import ProteomicsLibrary.DbTool;
import org.apache.commons.math3.util.Pair;
import proteomics.Index.BuildIndex;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

public class BuildDecoyProts implements Callable<BuildDecoyProts.Entry> {
    private final Map<String, String> parameterMap;

    private final BuildIndex buildIndex;

    private final String protId;

    public BuildDecoyProts(Map<String, String> parameterMap, BuildIndex buildIndex, String protId) throws Exception {
        this.parameterMap = parameterMap;
        this.buildIndex = buildIndex;
        this.protId = protId;
    }

    @Override
    public Entry call() throws Exception {
        boolean addDecoy = parameterMap.get("add_decoy").contentEquals("1");

        Entry entry = new Entry();
        entry.protId = protId;
        String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
        if (addDecoy) {
            String decoyProtSeq = DbTool.shuffleProtKeepKR(protSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1).replace('I', 'L');
            entry.decoyProtSeq = decoyProtSeq;
        }
        return entry;
    }

    public class Entry {
        String protId = "";
        String decoyProtSeq = "";

        Entry() {
        }
    }
}
