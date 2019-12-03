package FACSWebsiteEnd.Entity;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:58
 * QQ:776748935
 */
public class FileInfo {

    private String filename;
    private String path;
    private String fullpath;
    private String extension;

    public FileInfo() {
    }

    public FileInfo(String filename, String path, String fullpath, String extension) {
        this.filename = filename;
        this.path = path;
        this.fullpath = fullpath;
        this.extension = extension;
    }

    public String getFilename() {
        return filename;
    }

    public void setFilename(String filename) {
        this.filename = filename;
    }

    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getFullpath() {
        return fullpath;
    }

    public void setFullpath(String fullpath) {
        this.fullpath = fullpath;
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
                "filename='" + filename + '\'' +
                ", path='" + path + '\'' +
                ", fullpath='" + fullpath + '\'' +
                ", extension='" + extension + '\'' +
                '}';
    }
}
