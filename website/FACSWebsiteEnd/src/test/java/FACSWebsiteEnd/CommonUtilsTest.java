package FACSWebsiteEnd;

import FACSWebsiteEnd.utils.CommonUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

/**
 * @Author: HiramHe
 * @Date: 2019/12/8 16:13
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class CommonUtilsTest {

    @Test
    public void testGetCurrentTime(){
        System.out.println(CommonUtils.getCurrentTime());
    }


}
