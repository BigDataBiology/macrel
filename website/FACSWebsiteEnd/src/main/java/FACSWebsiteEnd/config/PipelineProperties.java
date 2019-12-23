package FACSWebsiteEnd.config;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.stereotype.Component;

/**
 * @Author: HiramHe
 * @Date: 2019/12/20 10:17
 * QQ:776748935
 */
@Component
@ConfigurationProperties(prefix = "pipeline")
public class PipelineProperties {
    //@Value("${pipeline.home}")
    private String home;
    //@Value("${pipeline.inputDir}")
    private String inputDir;
    //@Value("${pipeline.outputDir}")
    private String outputDir;

    public String getHome() {
        return home;
    }

    public void setHome(String home) {
        this.home = home;
    }

    public String getInputDir() {
        return inputDir;
    }

    public void setInputDir(String inputDir) {
        this.inputDir = inputDir;
    }

    public String getOutputDir() {
        return outputDir;
    }

    public void setOutputDir(String outputDir) {
        this.outputDir = outputDir;
    }

    @Override
    public String toString() {
        return "PipelineConfiguration{" +
                "home='" + home + '\'' +
                ", inputDir='" + inputDir + '\'' +
                ", outputDir='" + outputDir + '\'' +
                '}';
    }
}
