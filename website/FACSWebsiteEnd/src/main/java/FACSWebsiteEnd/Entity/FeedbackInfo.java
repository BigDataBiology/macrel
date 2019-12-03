package FACSWebsiteEnd.Entity;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 19:50
 * QQ:776748935
 */
public class FeedbackInfo {

    private Integer id;
    private String content;
    private String userEmail;

    public Integer getId() {
        return id;
    }

    public void setId(Integer id) {
        this.id = id;
    }

    public String getContent() {
        return content;
    }

    public void setContent(String content) {
        this.content = content;
    }

    public String getUserEmail() {
        return userEmail;
    }

    public void setUserEmail(String userEmail) {
        this.userEmail = userEmail;
    }

    @Override
    public String toString() {
        return "FeedbackInfo{" +
                "id=" + id +
                ", content='" + content + '\'' +
                ", userEmail='" + userEmail + '\'' +
                '}';
    }
}
