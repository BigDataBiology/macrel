package FACSWebsiteEnd.Entity;

import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/12/10 20:05
 * QQ:776748935
 */
public class PredictionOut {
    private String filePath;
    private List<Object> objects;

    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

    public List<Object> getObjects() {
        return objects;
    }

    public void setObjects(List<Object> objects) {
        this.objects = objects;
    }

    @Override
    public String toString() {
        return "PredictionOut{" +
                "filePath='" + filePath + '\'' +
                ", objects=" + objects +
                '}';
    }
}
