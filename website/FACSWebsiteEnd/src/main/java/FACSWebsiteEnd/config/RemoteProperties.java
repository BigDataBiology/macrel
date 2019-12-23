package FACSWebsiteEnd.config;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.stereotype.Component;

/**
 * @Author: HiramHe
 * @Date: 2019/12/19 21:14
 * QQ:776748935
 */

@Component
@ConfigurationProperties(prefix = "remote")
public class RemoteProperties {
    //@Value("${remote.enableRemote}")
    private Boolean enableRemote;
    //@Value("${remote.ip}")
    private String ip;
    //@Value("${remote.port}")
    private Integer port;
    //@Value("${remote.username}")
    private String username;
    //@Value("${remote.password}")
    private String password;

    public Boolean getEnableRemote() {
        return enableRemote;
    }

    public void setEnableRemote(Boolean enableRemote) {
        this.enableRemote = enableRemote;
    }

    public String getIp() {
        return ip;
    }

    public void setIp(String ip) {
        this.ip = ip;
    }

    public Integer getPort() {
        return port;
    }

    public void setPort(Integer port) {
        this.port = port;
    }

    public String getUsername() {
        return username;
    }

    public void setUsername(String username) {
        this.username = username;
    }

    public String getPassword() {
        return password;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    @Override
    public String toString() {
        return "RemoteConfiguration{" +
                "enableRemote=" + enableRemote +
                ", ip='" + ip + '\'' +
                ", port=" + port +
                ", username='" + username + '\'' +
                ", password='" + password + '\'' +
                '}';
    }
}
