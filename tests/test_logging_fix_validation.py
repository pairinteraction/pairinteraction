#!/usr/bin/env python3
"""
Integration test for spdlog logging fix (Issue #293).

This test reproduces the exact scenario described in the issue
and validates that the fix works as expected.
"""

import logging
import sys
from pathlib import Path


def test_logging_issue_fix():
    """
    Test that reproduces the issue from #293 and validates the fix.
    
    Issue: When creating Database instances, the log messages from C++
    don't appear in Python logging until another pairinteraction function
    is called.
    
    Root cause: _flush_pending_logs() was not being called after Database.__init__
    because the decorator only decorated methods that don't start with "__".
    
    Fix: Modified decorate_module_with_flush_logs() to also decorate __init__ methods.
    
    Expected behavior after fix: Log messages should appear immediately
    after Database constructor completes due to automatic log flushing.
    """
    # Set up Python logging as in the issue example
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    
    # Capture log messages for verification
    captured_logs = []
    
    class TestHandler(logging.Handler):
        def emit(self, record):
            captured_logs.append(record.getMessage())
    
    test_handler = TestHandler()
    logger.addHandler(test_handler)
    
    try:
        # Try to import pairinteraction
        import pairinteraction.real as pi
        
        print("Testing Database creation logging...")
        
        # Clear any existing logs
        captured_logs.clear()
        
        # Create databases - this should trigger log messages AND automatic flush
        # because Database.__init__ is now decorated with _flush_logs_after()
        db1 = pi.Database(database_dir="test1")
        
        # Check if logs appeared immediately after Database constructor
        database_logs_after_db1 = [log for log in captured_logs if "database" in log.lower()]
        
        db2 = pi.Database(database_dir="test2")
        
        database_logs_after_db2 = [log for log in captured_logs if "database" in log.lower()]
        
        if len(database_logs_after_db2) >= 2:
            print("✓ SUCCESS: Database creation logs appeared immediately")
            print(f"  Found {len(database_logs_after_db2)} database-related log messages")
            for log in database_logs_after_db2:
                print(f"  - {log}")
            return True
        elif len(database_logs_after_db1) >= 1:
            print("✓ PARTIAL SUCCESS: Some database logs appeared")
            print(f"  Found {len(database_logs_after_db1)} database-related log messages after first Database")
            print(f"  Found {len(database_logs_after_db2)} database-related log messages after second Database")
            return True
        else:
            print("✗ ISSUE: Database logs did not appear immediately")
            print(f"  Expected: 1+ database logs, Found: {len(database_logs_after_db2)}")
            
            # Try the workaround mentioned in the issue
            print("  Trying workaround: creating KetAtom to trigger flush...")
            pi.KetAtom("Rb", n=60, l=0, m=0.5)
            
            database_logs_after_workaround = [log for log in captured_logs if "database" in log.lower()]
            if len(database_logs_after_workaround) > len(database_logs_after_db2):
                print("  ⚠ Database logs appeared only after additional function call")
                print("  This indicates the issue is NOT fixed")
                return False
            else:
                print("  ✗ Database logs still didn't appear - different issue?")
                return False
            
    except ImportError as e:
        print(f"Cannot test - pairinteraction not available: {e}")
        print("This test requires pairinteraction to be built and installed")
        return None
    
    finally:
        logger.removeHandler(test_handler)


def main():
    print("=" * 70)
    print("INTEGRATION TEST: spdlog logging fix for Issue #293")
    print("Fix: Decorator now includes __init__ methods for log flushing")
    print("=" * 70)
    print()
    
    result = test_logging_issue_fix()
    
    print()
    print("=" * 70)
    if result is True:
        print("✓ TEST PASSED: Logging fix works correctly")
        exit_code = 0
    elif result is False:
        print("✗ TEST FAILED: Logging issue persists")
        exit_code = 1
    else:
        print("? TEST SKIPPED: Cannot test without pairinteraction")
        exit_code = 2
    print("=" * 70)
    
    return exit_code


if __name__ == "__main__":
    sys.exit(main())