package FACSWebsiteEnd.Entity;

import org.springframework.web.multipart.MultipartFile;

/**
 * @Author: HiramHe
 * @Date: 2019/11/29 11:18
 * QQ:776748935
 */
public class DataUploaded {

    private String sequence;
    private MultipartFile file;

    private String item1;
    private String item2;

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public MultipartFile getFile() {
        return file;
    }

    public void setFile(MultipartFile file) {
        this.file = file;
    }

    public String getItem1() {
        return item1;
    }

    public void setItem1(String item1) {
        this.item1 = item1;
    }

    public String getItem2() {
        return item2;
    }

    public void setItem2(String item2) {
        this.item2 = item2;
    }

    @Override
    public String toString() {
        return "DataUploaded{" +
                "sequence='" + sequence + '\'' +
                ", file=" + file +
                ", item1='" + item1 + '\'' +
                ", item2='" + item2 + '\'' +
                '}';
    }
}
