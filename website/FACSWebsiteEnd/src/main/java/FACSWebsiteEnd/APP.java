package FACSWebsiteEnd;

import FACSWebsiteEnd.config.PipelineProperties;
import FACSWebsiteEnd.config.RemoteProperties;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.EnableAutoConfiguration;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.context.properties.EnableConfigurationProperties;

import java.util.Scanner;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:56
 * QQ:776748935
 */

@SpringBootApplication
public class APP {

    public static void main(String[] args) {
        SpringApplication springApplication = new SpringApplication(APP.class);
        // 将命令行参数的优先级设置为最高。
        springApplication.setAddCommandLineProperties(true);
        springApplication.run(args);
    }
}
