package FACSWebsiteEnd.Entity;

/**
 * @Author: HiramHe
 * @Date: 2019/12/7 17:46
 * QQ:776748935
 */
public class FacsOutTsv {

    private String Access;
    private String Sequence;
    private String AMP_family;
    private Double AMP_probability;
    private String Hemolytic;
    private Double Hemolytic_probability;

    public FacsOutTsv() {
    }

    public FacsOutTsv(String access, String sequence, String AMP_family, Double AMP_probability, String hemolytic, Double hemolytic_probability) {
        Access = access;
        Sequence = sequence;
        this.AMP_family = AMP_family;
        this.AMP_probability = AMP_probability;
        Hemolytic = hemolytic;
        Hemolytic_probability = hemolytic_probability;
    }

    public String getAccess() {
        return Access;
    }

    public void setAccess(String access) {
        Access = access;
    }

    public String getSequence() {
        return Sequence;
    }

    public void setSequence(String sequence) {
        Sequence = sequence;
    }

    public String getAMP_family() {
        return AMP_family;
    }

    public void setAMP_family(String AMP_family) {
        this.AMP_family = AMP_family;
    }

    public Double getAMP_probability() {
        return AMP_probability;
    }

    public void setAMP_probability(Double AMP_probability) {
        this.AMP_probability = AMP_probability;
    }

    public String getHemolytic() {
        return Hemolytic;
    }

    public void setHemolytic(String hemolytic) {
        Hemolytic = hemolytic;
    }

    public Double getHemolytic_probability() {
        return Hemolytic_probability;
    }

    public void setHemolytic_probability(Double hemolytic_probability) {
        Hemolytic_probability = hemolytic_probability;
    }

    @Override
    public String toString() {
        return "FacsOutTsv{" +
                ", Access='" + Access + '\'' +
                ", Sequence='" + Sequence + '\'' +
                ", AMP_family='" + AMP_family + '\'' +
                ", AMP_probability=" + AMP_probability +
                ", Hemolytic='" + Hemolytic + '\'' +
                ", Hemolytic_probability=" + Hemolytic_probability +
                '}';
    }
}
