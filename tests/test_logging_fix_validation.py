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
    
    Expected behavior after fix: Log messages should appear immediately
    when flushed, not requiring additional function calls.
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
        
        # Create databases - this should trigger log messages
        db1 = pi.Database(database_dir="test1")
        db2 = pi.Database(database_dir="test2")
        
        # Check if logs appeared without needing additional function calls
        database_logs = [log for log in captured_logs if "database directory" in log.lower()]
        
        if len(database_logs) >= 2:
            print("✓ SUCCESS: Database creation logs appeared immediately")
            print(f"  Found {len(database_logs)} database-related log messages")
            for log in database_logs:
                print(f"  - {log}")
            return True
        else:
            print("✗ ISSUE: Database logs did not appear immediately")
            print(f"  Expected: 2+ database logs, Found: {len(database_logs)}")
            
            # Try the workaround mentioned in the issue
            print("  Trying workaround: creating KetAtom to trigger flush...")
            pi.KetAtom("Rb", n=60, l=0, m=0.5)
            
            database_logs_after = [log for log in captured_logs if "database directory" in log.lower()]
            if len(database_logs_after) > len(database_logs):
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