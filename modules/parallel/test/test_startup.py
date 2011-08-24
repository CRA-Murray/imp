import IMP
import IMP.test
import IMP.parallel
import util
import tasks

class StartupTests(IMP.test.TestCase):
    """Test startup of parallel jobs"""

    def test_simple(self):
        """Test that slave tasks can start up and communicate"""
        m = IMP.parallel.Manager()
        m.add_slave(IMP.parallel.LocalSlave())
        c = m.get_context()
        for i in range(10):
            c.add_task(tasks.SimpleTask(i))
        results = list(c.get_results_unordered())
        results.sort()
        self.assertEqual(results, list(range(10)))
        util.unlink("slave0.output")

    def test_startup_heartbeat(self):
        """Make sure that startup failures cause a timeout"""
        def empty_task():
            pass
        m = IMP.parallel.Manager(python="/path/does/not/exist")
        m.heartbeat_timeout = 0.1
        m.add_slave(IMP.parallel.LocalSlave())
        c = m.get_context()
        c.add_task(empty_task)
        self.assertRaises(IMP.parallel.NetworkError, list,
                          c.get_results_unordered())
        util.unlink('slave0.output')

if __name__ == '__main__':
    IMP.test.main()
