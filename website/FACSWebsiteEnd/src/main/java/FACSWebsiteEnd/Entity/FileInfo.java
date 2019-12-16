package FACSWebsiteEnd.Entity;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:58
 * QQ:776748935
 */
public class FileInfo {

    private String filenameWithExtension;
    private String filenameWithOutExtension;
    private String dir;
    private String path;
    private String extension;

    public FileInfo() {
    }

    public FileInfo(String filenameWithExtension, String filenameWithOutExtension, String dir, String path, String extension) {
        this.filenameWithExtension = filenameWithExtension;
        this.filenameWithOutExtension = filenameWithOutExtension;
        this.dir = dir;
        this.path = path;
        this.extension = extension;
    }

    public String getFilenameWithExtension() {
        return filenameWithExtension;
    }

    public void setFilenameWithExtension(String filenameWithExtension) {
        this.filenameWithExtension = filenameWithExtension;
    }

    public String getFilenameWithOutExtension() {
        return filenameWithOutExtension;
    }

    public void setFilenameWithOutExtension(String filenameWithOutExtension) {
        this.filenameWithOutExtension = filenameWithOutExtension;
    }

    public String getDir() {
        return dir;
    }

    public void setDir(String dir) {
        this.dir = dir;
    }

    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getExtension() {
        return extension;
    }

    public void setExtension(String extension) {
        this.extension = extension;
    }

    @Override
    public String toString() {
        return "FileInfo{" +
                "filenameWithExtension='" + filenameWithExtension + '\'' +
                ", filenameWithOutExtension='" + filenameWithOutExtension + '\'' +
                ", dir='" + dir + '\'' +
                ", path='" + path + '\'' +
                ", extension='" + extension + '\'' +
                '}';
    }
}
