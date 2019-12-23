package FACSWebsiteEnd.Entity;

import org.springframework.web.multipart.MultipartFile;

/**
 * @Author: HiramHe
 * @Date: 2019/11/29 11:18
 * QQ:776748935
 */
public class PredictionForm {

    private String textData;
    private MultipartFile file;

    private String dataType;

    public String getTextData() {
        return textData;
    }

    public void setTextData(String textData) {
        this.textData = textData;
    }

    public MultipartFile getFile() {
        return file;
    }

    public void setFile(MultipartFile file) {
        this.file = file;
    }

    public String getDataType() {
        return dataType;
    }

    public void setDataType(String dataType) {
        this.dataType = dataType;
    }

    @Override
    public String toString() {
        return "PredictionForm{" +
                "textData='" + textData + '\'' +
                ", file=" + file +
                ", dataType='" + dataType + '\'' +
                '}';
    }
}
