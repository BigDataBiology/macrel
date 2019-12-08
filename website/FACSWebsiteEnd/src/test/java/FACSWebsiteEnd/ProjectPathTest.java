package FACSWebsiteEnd;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.ClassUtils;

/**
 * @Author: HiramHe
 * @Date: 2019/12/7 14:06
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class ProjectPathTest {

    @Test
    public void testGetProjectAbsolutePath(){
        // 获取classes目录绝对路径
        String path = ClassUtils.getDefaultClassLoader().getResource("").getPath();
        System.out.println(path);

        // 获取当前项目路径的地址
        System.out.println(System.getProperty("user.dir"));
    }

}
