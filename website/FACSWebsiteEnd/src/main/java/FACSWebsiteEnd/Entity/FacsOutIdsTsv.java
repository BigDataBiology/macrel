package FACSWebsiteEnd.Entity;

/**
 * @Author: HiramHe
 * @Date: 2019/12/7 18:52
 * QQ:776748935
 */
public class FacsOutIdsTsv {
    private String smORF_num;
    private String Sequence;
    private String UniqueMark;

    public FacsOutIdsTsv() {
    }

    public FacsOutIdsTsv(String smORF_num, String sequence, String uniqueMark) {
        this.smORF_num = smORF_num;
        Sequence = sequence;
        UniqueMark = uniqueMark;
    }

    public String getSmORF_num() {
        return smORF_num;
    }

    public void setSmORF_num(String smORF_num) {
        this.smORF_num = smORF_num;
    }

    public String getSequence() {
        return Sequence;
    }

    public void setSequence(String sequence) {
        Sequence = sequence;
    }

    public String getUniqueMark() {
        return UniqueMark;
    }

    public void setUniqueMark(String uniqueMark) {
        UniqueMark = uniqueMark;
    }

    @Override
    public String toString() {
        return "FacsOutIdsTsv{" +
                "smORF_num='" + smORF_num + '\'' +
                ", Sequence='" + Sequence + '\'' +
                ", UniqueMark='" + UniqueMark + '\'' +
                '}';
    }
}
