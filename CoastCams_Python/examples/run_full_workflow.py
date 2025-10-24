#!/usr/bin/env python3
"""
Example: Run Full CoastCams Workflow

This script demonstrates how to run the complete CoastCams analysis
workflow with default or custom settings.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coastcams.config import CoastCamsConfig
from main import CoastCamsWorkflow


def main():
    """Run full workflow example."""

    print("="*70)
    print("CoastCams - Full Workflow Example")
    print("="*70)

    # Option 1: Run with config.yaml in parent directory
    print("\nOption 1: Running with config.yaml...")
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'config.yaml')
    workflow = CoastCamsWorkflow(config_file=config_path)
    workflow.run_full_analysis()

    # Option 2: Run with custom configuration
    # Uncomment to use custom config file
    # workflow = CoastCamsWorkflow(config_file='my_config.yaml')
    # workflow.run_full_analysis()

    # Option 3: Customize settings programmatically
    # config = CoastCamsConfig()
    # config.shoreline_method = 2  # Use red-minus-blue method
    # config.wave_period_min = 5.0
    # config.wave_period_max = 20.0
    # # Save customized config
    # config.save_to_file('custom_config.yaml')


if __name__ == '__main__':
    main()
